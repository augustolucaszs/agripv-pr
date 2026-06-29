## Performance Ratio (PR) for AGRIPV — Florianópolis, Brazil
#
# OUTPUT TABLES (in fotovoltaica_clean DB):
#   AGRIPV_raw_pr_mppt_01..07  per-minute: TIMESTAMP, mppt_id, power_w, gpoa, ghi, position_deg, ...
#   AGRIPV_daily_pr            daily IEC 61724 PR: date, mppt_id, pr, energy_kwh, poa_kwh_m2, ...
#
# IEC 61724 PR formula:
#   PR = Σ(power_W) / (Σ(GPOA) / 1000 × P_nom_W)   — Δt cancels (constant 1-min interval)
#   Included minutes: GHI ≥ 50 W/m², power_w > 0, GPOA > 0, and not tracker_lag_flag.
#   Computed twice per MPPT per day: full day, and 10:00-14:00 local-time peak window
#   (pr_peak / energy_kwh_peak / poa_kwh_m2_peak / n_minutes_peak).


# bring power and irrad gpoa ghi to PA dash. add selection item possibility
# add system status for HOYMILES
# add power/m2

# nth
# PR temperature corrected
# Bifacial PR

import argparse
import logging
import os
from pathlib import Path
from datetime import datetime, timedelta

from dotenv import load_dotenv
from sqlalchemy import create_engine, text
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
    1: dict(nominal_w=10695, tracking=False, tilt=27.0, azimuth=0.0),
    2: dict(nominal_w=9350,  tracking=True),
    3: dict(nominal_w=8800,  tracking=True),
    4: dict(nominal_w=8800,  tracking=True),
    5: dict(nominal_w=9350,  tracking=True),
    6: dict(nominal_w=12250, tracking=True),
    7: dict(nominal_w=12250, tracking=True),
}

# DC string(s) per MPPT (1-based → dcw_N columns)
MPPT_STRINGS = {1: [1], 2: [2], 3: [3], 4: [4], 5: [5], 6: [6], 7: [7]}

# Tracker groups (trackers_suntrack_tcu.group_id) per MPPT
# Mapping adjusted for validation/comparison as provided by field notes.
MPPT_GROUPS = {2: [9], 3: [8], 4: [7], 5: [6], 6: [4, 5], 7: [2, 3]}

GHI_MIN = 50  # W/m² — (legacy) exclude low-irradiance minutes
GPOA_MIN = 50  # W/m² — single PR filter: only minutes with GPOA above this are kept

TRACKER_SMOOTH_WIN  = 3    # minutes, centered rolling mean applied to tracker angle for GPOA
TRACKER_LAG_MAX_DEG = 2.0  # |position - target| above this = tracker in transit, unreliable minute

PEAK_HOUR_START, PEAK_HOUR_END = 10, 14  # local time [10:00, 14:00) — secondary "peak window" PR


def _fetch(engine, table, columns, start, end):
    cols = ', '.join(f'"{c}"' for c in columns)
    q = (f'SELECT "TIMESTAMP", {cols} FROM "{table}" '
         f"WHERE \"TIMESTAMP\" BETWEEN '{start}' AND '{end}' "
         f"ORDER BY \"TIMESTAMP\"")
    with engine.connect() as conn:
        return pd.read_sql(q, conn)


def _compute_poa(df, tilt, azimuth_series):
    """
    Compute plane-of-array irradiance via Erbs decomposition + Hay-Davies model.
    df must be DatetimeIndex with a 'ghi' column.
    tilt: scalar (degrees) or column name string → abs() applied.
    Returns ghi_clipped, poa_global (both Series aligned to df.index).
    """
    ghi = df['ghi'].clip(lower=0)
    low = ghi < 50
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
        model='haydavies',
    )
    return ghi_c, poa['poa_global']


def _pr_stats(filt, nominal_w):
    """IEC 61724 PR + energy/POA sums for an already-filtered per-minute slice."""
    n = len(filt)
    if n == 0:
        return dict(pr=None, energy_kwh=0.0, poa_kwh_m2=0.0, n_minutes=0)

    # 1-min data: Σ(W) / 60 = Wh; /1000 = kWh
    energy_kwh = float(filt['power_w'].sum()) / 60 / 1000
    poa_kwh_m2 = float(filt['gpoa'].sum())    / 60 / 1000
    pr = energy_kwh / (poa_kwh_m2 * nominal_w / 1000) if poa_kwh_m2 > 0 else None

    return dict(
        pr=round(pr, 4) if pr is not None else None,
        energy_kwh=round(energy_kwh, 4),
        poa_kwh_m2=round(poa_kwh_m2, 4),
        n_minutes=n,
    )


def process_day(start_date, end_date):
    """
    Returns
    -------
    mppt_dfs : dict[int → DataFrame]  per-minute measurements per MPPT
    daily_df : DataFrame              one row per MPPT — IEC 61724 daily PR
    """
    logging.info(f"Processing {start_date}")

    inv_raw = _fetch(engine_source, 'agripv_inverter_sma',
                     [f'dcw_{i}' for i in range(1, 8)], start_date, end_date)
    met_raw = _fetch(engine_source, 'SapAlbedo_1m',
                     ['GHIA_SMP22_Comp_Avg'], start_date, end_date)
    trk_raw = _fetch(engine_source, 'trackers_suntrack_tcu',
                     ['position_a1_degree', 'targetangle_a1_degree', 'btactive_a1', 'group_id'],
                     start_date, end_date)

    inv_raw['TIMESTAMP'] = pd.to_datetime(inv_raw['TIMESTAMP'])
    met_raw['TIMESTAMP'] = pd.to_datetime(met_raw['TIMESTAMP'])
    trk_raw['TIMESTAMP'] = pd.to_datetime(trk_raw['TIMESTAMP'])

    inv = inv_raw.set_index('TIMESTAMP').sort_index()
    for i in range(1, 8):
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
        df1['position_deg']    = MPPT[1]['tilt']
        df1['target_deg']      = MPPT[1]['tilt']
        df1['tracker_active']  = False
        df1['tracker_lag_flag']= False

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

        # Tracker lags its setpoint by a few degrees for ~1 min after leaving stow or
        # during fast backtracking corrections — the angle reading for that minute is
        # "in transit", not a stable orientation. Flag it so it's excluded from PR sums.
        df['tracker_lag_flag'] = (df['position_deg'] - df['target_deg']).abs() > TRACKER_LAG_MAX_DEG

        # Smooth the raw angle before computing GPOA to reduce minute-to-minute noise.
        # The smoothed angle is only used for the GPOA calc — stored position_deg stays raw.
        df['position_smooth'] = df['position_deg'].rolling(
            TRACKER_SMOOTH_WIN, center=True, min_periods=1).mean()

        # N-S HSAT azimuth: East-facing (negative angle) = 90°, West-facing (positive) = 270°.
        # Using solar azimuth here would treat the tracker as 2-axis, overestimating GPOA.
        hsat_az = pd.Series(
            np.where(df['position_smooth'] >= 0, 270.0, 90.0), index=df.index
        )
        ghi_c, gpoa = _compute_poa(df, tilt='position_smooth', azimuth_series=hsat_az)
        df['ghi']  = ghi_c
        df['gpoa'] = gpoa
        df.drop(columns=['position_smooth'], inplace=True)

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
        mask = (
            (df['gpoa'].fillna(0) > GPOA_MIN)
            & df['power_w'].notna() & (df['power_w'] > 0)
            & ~df['tracker_lag_flag']
        )

        local_hour = df['TIMESTAMP'].dt.tz_localize('UTC').dt.tz_convert(SITE.tz).dt.hour
        peak_mask = mask & local_hour.between(PEAK_HOUR_START, PEAK_HOUR_END - 1)

        full = _pr_stats(df[mask], nominal_w)
        peak = _pr_stats(df[peak_mask], nominal_w)

        daily_rows.append(dict(
            date=start_date, mppt_id=mid,
            pr=full['pr'], energy_kwh=full['energy_kwh'], poa_kwh_m2=full['poa_kwh_m2'],
            n_minutes=full['n_minutes'],
            pr_peak=peak['pr'], energy_kwh_peak=peak['energy_kwh'],
            poa_kwh_m2_peak=peak['poa_kwh_m2'], n_minutes_peak=peak['n_minutes'],
            nominal_kw=nominal_w / 1000,
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
    parser.add_argument('start_date', nargs='?', default=None,
                        help='YYYY-MM-DD (default: yesterday, for cron use)')
    parser.add_argument('end_date',   nargs='?', default=None,
                        help='YYYY-MM-DD (default: same as start_date)')
    parser.add_argument('--csv',        action='store_true',
                        help='Write CSV files instead of pushing to DB')
    parser.add_argument('--daily-only', action='store_true',
                        help='Push only the daily PR table; skip per-minute tables')
    parser.add_argument('-d', '--debug',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO')
    args = parser.parse_args()

    if args.start_date is None:
        # No args (e.g. cron at midnight) → process yesterday, local time.
        yesterday = (pd.Timestamp.now(tz=SITE.tz) - pd.Timedelta(days=1)).strftime('%Y-%m-%d')
        args.start_date = yesterday
    if args.end_date is None:
        args.end_date = args.start_date

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
            # Idempotent rerun: clear the day first so the TIMESTAMP primary key
            # on AGRIPV_raw_pr_mppt_* never raises a duplicate-key error.
            with engine_clean.begin() as conn:
                if not args.daily_only:
                    for mid in mppt_dfs:
                        conn.execute(text(
                            f'DELETE FROM "AGRIPV_raw_pr_mppt_{mid:02d}" '
                            'WHERE "TIMESTAMP" >= :s AND "TIMESTAMP" < :e'),
                            {'s': ds, 'e': de})
                conn.execute(text(
                    'DELETE FROM "AGRIPV_daily_pr" WHERE date = :d'), {'d': ds})

            if not args.daily_only:
                for mid, df in mppt_dfs.items():
                    df.to_sql(f'AGRIPV_raw_pr_mppt_{mid:02d}', engine_clean,
                              if_exists='append', index=False)
            daily_df.to_sql('AGRIPV_daily_pr', engine_clean,
                            if_exists='append', index=False)
            logging.info(f"  DB push done for {ds}")


if __name__ == '__main__':
    main()
