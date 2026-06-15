## Performance Ratio (PR) for AGRIPV — Florianópolis, Brazil
#
# OUTPUT TABLES (in fotovoltaica_clean DB):
#   AGRIPV_raw_pr_mppt_01..07  per-minute: TIMESTAMP, mppt_id, power_w, gpoa, ghi, position_deg, ...
#   AGRIPV_daily_pr            daily IEC 61724 PR: date, mppt_id, pr, energy_kwh, poa_kwh_m2, ...
#
# IEC 61724 PR formula:
#   PR = Σ(power_W) / (Σ(GPOA) / 1000 × P_nom_W)   — Δt cancels (constant 1-min interval)
#   Only minutes with GHI ≥ 50 W/m² are included.

import argparse
import logging
import os
from pathlib import Path
from datetime import datetime, timedelta

from dotenv import load_dotenv
from sqlalchemy import create_engine
import numpy as np
import pandas as pd
from pvlib import location, irradiance

load_dotenv()

_src_url = os.environ.get("DB_SOURCE_URL")
_cln_url = os.environ.get("DB_CLEAN_URL")
if not _src_url or not _cln_url:
    raise EnvironmentError("DB_SOURCE_URL and DB_CLEAN_URL must be set in .env")

engine_source = create_engine(_src_url, echo=False)
engine_clean  = create_engine(_cln_url, echo=False)

SITE = location.Location(-27.597, -48.549, tz='America/Sao_Paulo', altitude=3)

MPPT = {
    1: dict(nominal_w=10695, tracking=False, tilt=15.0, azimuth=0.0),
    2: dict(nominal_w=9265,  tracking=True),
    3: dict(nominal_w=8720,  tracking=True),
    4: dict(nominal_w=8720,  tracking=True),
    5: dict(nominal_w=9265,  tracking=True),
    6: dict(nominal_w=12250, tracking=True),
    7: dict(nominal_w=12250, tracking=True),
}

# DC string(s) per MPPT (1-based → dcw_N columns)
MPPT_STRINGS = {1: [1], 2: [2], 3: [3], 4: [4], 5: [5], 6: [6], 7: [7]}

# Tracker groups (trackers_suntrack_tcu.group_id) per MPPT
# MPPTs 6 & 7 each have 2 groups — average their positions for a single GPOA
MPPT_GROUPS = {2: [2], 3: [3], 4: [4], 5: [5], 6: [6, 7], 7: [8, 9]}

GHI_MIN = 50  # W/m² — exclude low-irradiance minutes from PR


def _fetch(engine, table, columns, start, end):
    cols = ', '.join(f'"{c}"' for c in columns)
    q = (f'SELECT "TIMESTAMP", {cols} FROM "{table}" '
         f"WHERE \"TIMESTAMP\" BETWEEN '{start}' AND '{end}' "
         f"ORDER BY \"TIMESTAMP\"")
    with engine.connect() as conn:
        return pd.read_sql(q, conn)


def _compute_poa(df, tilt, azimuth_series):
    """
    Compute plane-of-array irradiance via Erbs decomposition + Perez model.
    df must be DatetimeIndex with a 'ghi' column.
    tilt: scalar (degrees) or column name string → abs() applied.
    Returns ghi_clipped, poa_global (both Series aligned to df.index).
    """
    ghi = df['ghi'].clip(lower=0)
    low = ghi < 20
    ghi_c = ghi.where(~low, 0.0)

    solar = SITE.get_solarposition(df.index)
    erbs  = irradiance.erbs(ghi=ghi_c, zenith=solar['zenith'], datetime_or_doy=df.index)
    dni   = erbs['dni'].where(~low, 0.0)
    dhi   = erbs['dhi'].where(~low, 0.0)

    surf_tilt = float(tilt) if isinstance(tilt, (int, float)) else df[tilt].astype(float).abs()

    poa = irradiance.get_total_irradiance(
        surface_tilt    = surf_tilt,
        surface_azimuth = azimuth_series,
        solar_zenith    = solar['zenith'],
        solar_azimuth   = solar['azimuth'],
        dni=dni, ghi=ghi_c, dhi=dhi,
        dni_extra=irradiance.get_extra_radiation(df.index),
        model='perez',
    )
    return ghi_c, poa['poa_global']


def process_day(start_date, end_date):
    """
    Returns
    -------
    mppt_dfs : dict[int → DataFrame]  per-minute measurements per MPPT
    daily_df : DataFrame              one row per MPPT — IEC 61724 daily PR
    """
    logging.info(f"Processing {start_date}")

    inv_raw = _fetch(engine_source, 'agripv_inverter_sma',
                     [f'dcw_{i}' for i in range(1, 13)], start_date, end_date)
    met_raw = _fetch(engine_source, 'SapAlbedo_1m',
                     ['GHIA_SMP22_Comp_Avg'], start_date, end_date)
    trk_raw = _fetch(engine_source, 'trackers_suntrack_tcu',
                     ['position_a1_degree', 'targetangle_a1_degree', 'btactive_a1', 'group_id'],
                     start_date, end_date)

    inv_raw['TIMESTAMP'] = pd.to_datetime(inv_raw['TIMESTAMP'])
    met_raw['TIMESTAMP'] = pd.to_datetime(met_raw['TIMESTAMP'])
    trk_raw['TIMESTAMP'] = pd.to_datetime(trk_raw['TIMESTAMP'])

    inv = inv_raw.set_index('TIMESTAMP').sort_index()
    for i in range(1, 13):
        inv[f'dcw_{i}'] = pd.to_numeric(inv[f'dcw_{i}'], errors='coerce')

    met = met_raw.set_index('TIMESTAMP').sort_index()
    met['ghi'] = pd.to_numeric(met['GHIA_SMP22_Comp_Avg'], errors='coerce').clip(lower=0)

    for col in ['position_a1_degree', 'targetangle_a1_degree']:
        trk_raw[col] = pd.to_numeric(trk_raw[col], errors='coerce')
    trk_raw['btactive_a1'] = pd.to_numeric(trk_raw['btactive_a1'], errors='coerce').astype(bool)

    mppt_dfs: dict = {}

    # ── Fixed MPPT 1 (north-facing, 15° tilt) ────────────────────────────────
    df1 = met[['ghi']].copy()
    if not df1.empty:
        solar1 = SITE.get_solarposition(df1.index)
        ghi1, gpoa1 = _compute_poa(
            df1, tilt=MPPT[1]['tilt'],
            azimuth_series=pd.Series(MPPT[1]['azimuth'], index=df1.index),
        )
        df1['ghi']           = ghi1
        df1['gpoa']          = gpoa1
        df1['position_deg']  = MPPT[1]['tilt']
        df1['target_deg']    = MPPT[1]['tilt']
        df1['tracker_active']= False

        pwr = [f'dcw_{i}' for i in MPPT_STRINGS[1] if f'dcw_{i}' in inv.columns]
        if pwr:
            df1 = df1.join(inv[pwr], how='left')
            df1['power_w'] = df1[pwr].sum(axis=1)
            df1.drop(columns=pwr, inplace=True)
        else:
            df1['power_w'] = np.nan

        df1['mppt_id'] = 1
        mppt_dfs[1] = df1.reset_index()

    # ── Tracked MPPTs 2-7 ────────────────────────────────────────────────────
    for mid in range(2, 8):
        groups = MPPT_GROUPS[mid]
        grp_data = trk_raw[trk_raw['group_id'].isin(groups)].copy()
        if grp_data.empty:
            logging.warning(f"  No tracker data for MPPT{mid}")
            continue

        # For split MPPTs (6 & 7), average position across 2 groups → single GPOA
        grp_avg = (
            grp_data.groupby('TIMESTAMP')
            .agg(
                position_deg  =('position_a1_degree',    'mean'),
                target_deg    =('targetangle_a1_degree', 'mean'),
                tracker_active=('btactive_a1',           'any'),
            )
        )

        df = grp_avg.join(met[['ghi']], how='inner').sort_index()
        if df.empty:
            continue

        solar = SITE.get_solarposition(df.index)
        ghi_c, gpoa = _compute_poa(df, tilt='position_deg', azimuth_series=solar['azimuth'])
        df['ghi']  = ghi_c
        df['gpoa'] = gpoa

        pwr = [f'dcw_{i}' for i in MPPT_STRINGS[mid] if f'dcw_{i}' in inv.columns]
        if pwr:
            df = df.join(inv[pwr], how='left')
            df['power_w'] = df[pwr].sum(axis=1)
            df.drop(columns=pwr, inplace=True)
        else:
            df['power_w'] = np.nan

        df['mppt_id'] = mid
        mppt_dfs[mid] = df.reset_index()

    # ── IEC 61724 daily PR per MPPT ───────────────────────────────────────────
    daily_rows = []
    for mid, df in mppt_dfs.items():
        nominal_w = MPPT[mid]['nominal_w']
        mask = (df['ghi'] >= GHI_MIN) & df['power_w'].notna() & (df['gpoa'].fillna(0) > 0)
        filt = df[mask]
        n = len(filt)

        if n == 0:
            daily_rows.append(dict(date=start_date, mppt_id=mid, pr=None,
                                   energy_kwh=0.0, poa_kwh_m2=0.0,
                                   nominal_kw=nominal_w / 1000, n_minutes=0,
                                   tracking=MPPT[mid]['tracking']))
            continue

        # 1-min data: Σ(W) / 60 = Wh; /1000 = kWh
        energy_kwh  = float(filt['power_w'].sum()) / 60 / 1000
        poa_kwh_m2  = float(filt['gpoa'].sum())    / 60 / 1000
        pr = energy_kwh / (poa_kwh_m2 * nominal_w / 1000) if poa_kwh_m2 > 0 else None

        daily_rows.append(dict(
            date=start_date, mppt_id=mid,
            pr=round(pr, 4) if pr is not None else None,
            energy_kwh=round(energy_kwh, 4),
            poa_kwh_m2=round(poa_kwh_m2, 4),
            nominal_kw=nominal_w / 1000,
            n_minutes=n,
            tracking=MPPT[mid]['tracking'],
        ))

    daily_df = pd.DataFrame(daily_rows)
    logging.info("  PR: " + "  ".join(
        (f"MPPT{r['mppt_id']}={r['pr']:.3f}" if r['pr'] is not None
         else f"MPPT{r['mppt_id']}=N/A")
        for _, r in daily_df.iterrows()
    ))
    return mppt_dfs, daily_df


def main():
    parser = argparse.ArgumentParser(description='AGRIPV daily PR — IEC 61724')
    parser.add_argument('start_date', help='YYYY-MM-DD')
    parser.add_argument('end_date',   help='YYYY-MM-DD')
    parser.add_argument('--csv',        action='store_true',
                        help='Write CSV files instead of pushing to DB')
    parser.add_argument('--daily-only', action='store_true',
                        help='Push only the daily PR table; skip per-minute tables')
    parser.add_argument('-d', '--debug',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO')
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.debug),
                        format='%(asctime)s %(levelname)s %(message)s')

    dates = pd.date_range(
        start=datetime.strptime(args.start_date, '%Y-%m-%d'),
        end=  datetime.strptime(args.end_date,   '%Y-%m-%d'),
        freq='D',
    )

    for d in dates:
        ds = d.strftime('%Y-%m-%d')
        de = (d + timedelta(days=1)).strftime('%Y-%m-%d')
        mppt_dfs, daily_df = process_day(ds, de)

        if args.csv:
            out = Path(__file__).resolve().parent
            for mid, df in mppt_dfs.items():
                df.to_csv(out / f'mppt_{mid}_{ds}.csv', index=False)
            daily_df.to_csv(out / f'daily_pr_{ds}.csv', index=False)
            logging.info(f"  CSV written for {ds}")
        else:
            if not args.daily_only:
                for mid, df in mppt_dfs.items():
                    df.to_sql(f'AGRIPV_raw_pr_mppt_{mid:02d}', engine_clean,
                              if_exists='append', index=False)
            daily_df.to_sql('AGRIPV_daily_pr', engine_clean,
                            if_exists='append', index=False)
            logging.info(f"  DB push done for {ds}")


if __name__ == '__main__':
    main()
