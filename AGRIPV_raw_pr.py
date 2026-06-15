## FLOW OF THE SCRIPT - Performance Ratio (PR) calculation for PV devices
# 1. Connect to the PostgreSQL database using SQLAlchemy.
# 2. Define a function `get_inverter_data()` to retrieve necessary data for PR calculation.
# 3. Inside the function, execute SQL queries to fetch:
#    - Energy produced by each MPPT (from DC string power, in W).
#    - Plane-of-array irradiance (GPOA, in W/m²) computed from tracker angles + GHI.
# 4. Compute expected energy per square meter (kWh/m²) from GPOA.
#    (Area is not used in the current calculation; it can be applied later to get a true PR/efficiency.)
# 5. Compute a raw PR-like metric as energy produced divided by incident energy per m².
# 6. Handle any potential exceptions during database operations and ensure proper resource management.

import argparse
import logging
import os
from pathlib import Path
from datetime import datetime, timedelta

from dotenv import load_dotenv
from sqlalchemy import create_engine, inspect, text
import numpy as np
import pvlib
import pandas as pd
from pvlib import location, irradiance

load_dotenv()

_src_url = os.environ.get("DB_SOURCE_URL")
_cln_url = os.environ.get("DB_CLEAN_URL")
if not _src_url or not _cln_url:
    raise EnvironmentError(
        "DB_SOURCE_URL and DB_CLEAN_URL must be set in .env or environment"
    )

engine_source = create_engine(_src_url, echo=False)
engine_clean  = create_engine(_cln_url, echo=False)


metadata_dict = {
    'inverter_table': 'agripv_inverter_sma',
    # We need total AC power plus per-string DC power (dcw_n) to estimate per-MPPT output
    'inverter_columns': ['watts'] + [f'dcw_{i}' for i in range(1, 13)],
    'tracker_table': 'trackers_suntrack_tcu',
    'tracker_columns': ['position_a1_degree', 'targetangle_a1_degree', 'btactive_a1', "group_id"],
    'meteo_table': 'SapAlbedo_1m',
    'meteo_columns': ['GHIA_SMP22_Comp_Avg']
}

inverter_dict = {
    'mppts_power': [10695, 9265, 8720, 8720, 9265, 12250, 12250],
    'modules_per_mppt': [31, 17, 16, 16, 17, 35, 35],
    'module_power': [345, 545, 545, 545, 545, 350, 350],  # in W
    'modules_dummy': [1, 0, 1, 1, 0, 1, 1],
    'tracking_status': [False, True, True, True, True, True, True],
    'mppt_isc': [18.7, 14.91, 14.91, 14.91, 14.91, 18.75, 18.75],  # in A,
    'module_area': [2.173572, 2.59119, 2.59119, 2.59119, 2.59119, 2.173572, 2.173572],  # in m²
}

# Mapping of MPPT index (1-based) to the inverter string power columns that belong to it.
# This is used to sum DC power from strings associated with each MPPT for PR calculation.
# Adjust this mapping based on your string-to-MPPT wiring.
mppt_string_map = {
    1: [1],
    2: [2],
    3: [3],
    4: [4],
    5: [5],
    6: [6],
    7: [7],
    # strings 8 and 12 are left unused in this map; they are empty
}

def get_data_from_db(engine, table_name, columns, start_date, end_date):
    cols_sql = ', '.join(f'"{c}"' for c in columns)

    query = f"""
    SELECT
        "TIMESTAMP",
        {cols_sql}
    FROM
        "{table_name}"
    WHERE
        "TIMESTAMP" BETWEEN '{start_date}' AND '{end_date}'
    ORDER BY "TIMESTAMP" DESC
    """

    logging.debug(f"Executing query: {query}")

    with engine.connect() as conn:
        df = pd.read_sql(query, conn)

    return df


def process_day(start_date, end_date):
    logging.info(f"Processing data from {start_date} to {end_date}")

    inverter_data = get_data_from_db(engine_source, metadata_dict['inverter_table'], metadata_dict['inverter_columns'], start_date, end_date)
    logging.debug(f"Inverter data shape: {inverter_data.shape}")

    meteo_data = get_data_from_db(engine_source, metadata_dict['meteo_table'], metadata_dict['meteo_columns'], start_date, end_date)
    # Ensure numeric types
    meteo_data['GHIA_SMP22_Comp_Avg'] = pd.to_numeric(meteo_data['GHIA_SMP22_Comp_Avg'], errors='coerce')
    logging.debug(f"Meteo data shape: {meteo_data.shape}")

    tracker_data = get_data_from_db(engine_source, metadata_dict['tracker_table'], metadata_dict['tracker_columns'], start_date, end_date)
    # Ensure numeric types
    tracker_data['position_a1_degree'] = pd.to_numeric(tracker_data['position_a1_degree'], errors='coerce')
    tracker_data['targetangle_a1_degree'] = pd.to_numeric(tracker_data['targetangle_a1_degree'], errors='coerce')
    tracker_data['btactive_a1'] = pd.to_numeric(tracker_data['btactive_a1'], errors='coerce').astype(bool)
    # Only keep tracker groups that map to MPPTs (groups 2..9)
    tracker_data = tracker_data[tracker_data['group_id'].isin(range(2, 10))]
    logging.debug(f"Tracker data shape: {tracker_data.shape}")

    # Map group_id -> mppt_id (based on system wiring)
    group_to_mppt = {
        2: 2,
        3: 3,
        4: 4,
        5: 5,
        6: 6,
        7: 6,
        8: 7,
        9: 7,
    }

    # Build string-to-group mapping for split MPPTs
    # Note: adjust these assignments if your wiring differs.
    mppt_to_group = {
        2: [2],
        3: [3],
        4: [4],
        5: [5],
        6: [6, 7],
        7: [8, 9],
    }

    # Build group->string index map based on mppt_string_map and mppt_to_group
    group_string_map = {}
    for mppt_id, string_idx in mppt_string_map.items():
        groups = mppt_to_group.get(mppt_id, [])
        if not groups:
            continue
        if len(groups) == 1:
            group_string_map[groups[0]] = string_idx
        else:
            # split string indices evenly between groups
            half = len(string_idx) // len(groups)
            for i, grp in enumerate(groups):
                start = i * half
                end = (i + 1) * half if i < len(groups) - 1 else len(string_idx)
                group_string_map[grp] = string_idx[start:end]

    # Add the fixed MPPT (mppt 1) with assumed tilt 15 and north-facing
    fixed_mppt_id = 1
    fixed_tilt = 15.0
    fixed_azimuth = 0.0

    # Prepare inverter time series for power analysis
    inverter_data['TIMESTAMP'] = pd.to_datetime(inverter_data['TIMESTAMP'])
    inverter_data = inverter_data.set_index('TIMESTAMP').sort_index()

    # Prepare per-MPPT power from string DC power
    inverter_data = inverter_data.copy()
    for mppt_id, strings in mppt_string_map.items():
        cols = [f"dcw_{i}" for i in strings if f"dcw_{i}" in inverter_data.columns]
        if cols:
            inverter_data[f"mppt_{mppt_id}_power_w"] = inverter_data[cols].sum(axis=1)

    # Create per-group time series (tracker groups 2..9)
    per_group_dfs = {}
    for group_id in sorted(tracker_data['group_id'].unique()):
        group_df = tracker_data[tracker_data['group_id'] == group_id].copy()
        group_df['TIMESTAMP'] = pd.to_datetime(group_df['TIMESTAMP'])
        group_df = group_df.set_index('TIMESTAMP').sort_index()

        # Merge with meteo to compute irradiance
        merged = group_df.join(meteo_data.set_index(pd.to_datetime(meteo_data['TIMESTAMP'])), how='inner')
        merged = merged.sort_index()

        # Calculate solar position + POA (tracked)
        lat, lon = -27.597, -48.549
        tz = 'America/Sao_Paulo'
        altitude = 3
        site = location.Location(lat, lon, tz=tz, altitude=altitude, name='Florianopolis')
        solar_pos = site.get_solarposition(merged.index)
        dni_dhi = irradiance.erbs(
            ghi=merged['GHIA_SMP22_Comp_Avg'],
            zenith=solar_pos['zenith'],
            datetime_or_doy=merged.index,
        )
        merged['dni'] = dni_dhi['dni']
        merged['dhi'] = dni_dhi['dhi']
        merged.loc[merged['GHIA_SMP22_Comp_Avg'] < 20, 'GHIA_SMP22_Comp_Avg'] = 0
        merged['dni'] = merged['dni'].where(merged['GHIA_SMP22_Comp_Avg'] > 0, 0)
        merged['dhi'] = merged['dhi'].where(merged['GHIA_SMP22_Comp_Avg'] > 0, 0)

        poa = irradiance.get_total_irradiance(
            surface_tilt=merged['position_a1_degree'].abs(),
            surface_azimuth=solar_pos['azimuth'],
            solar_zenith=solar_pos['zenith'],
            solar_azimuth=solar_pos['azimuth'],
            dni=merged['dni'],
            ghi=merged['GHIA_SMP22_Comp_Avg'],
            dhi=merged['dhi'],
            dni_extra=irradiance.get_extra_radiation(merged.index),
            model='perez',
        )
        merged['gpoa'] = poa['poa_global']

        # Add mppt mapping
        mppt_id = group_to_mppt.get(int(group_id)) if pd.notna(group_id) else None
        merged['mppt_id'] = mppt_id
        merged['group_id'] = group_id

        # Compute group power from the associated strings (if any)
        string_idxs = group_string_map.get(int(group_id), [])
        cols = [f"dcw_{i}" for i in string_idxs if f"dcw_{i}" in inverter_data.columns]
        if cols:
            # Align inverter power to timestamps
            merged = merged.join(inverter_data[cols], how='left')
            merged[f'power_w'] = merged[cols].sum(axis=1)
        else:
            merged['power_w'] = None

        # Compute PR (power-based) using nominal MPPT rating
        if mppt_id is not None and mppt_id <= len(inverter_dict['modules_per_mppt']):
            module_count = inverter_dict['modules_per_mppt'][mppt_id - 1]
            module_area_m2 = inverter_dict['module_area'][mppt_id - 1]
            # If this mppt is split across two groups, assume each group is half the area
            area_m2 = module_count * module_area_m2
            if mppt_id in (6, 7):
                area_m2 = area_m2 / 2

            merged['area_m2'] = area_m2
            merged['mppt_nominal_w'] = inverter_dict['mppts_power'][mppt_id - 1]

            # expected power (W) = MPPT nominal power * (irradiance / 1000)
            # (irradiance is in W/m²; dividing by 1000 converts to per-unit of STC)
            irradiance_pu = merged['gpoa'] / 1000
            merged['expected_power_w'] = irradiance_pu * merged['mppt_nominal_w']

            mask = (merged['GHIA_SMP22_Comp_Avg'] >= 50) & (merged['expected_power_w'] > 0)
            merged['raw_pr'] = None
            merged.loc[mask, 'raw_pr'] = merged.loc[mask, 'power_w'] / merged.loc[mask, 'expected_power_w']
        else:
            merged['area_m2'] = None
            merged['mppt_nominal_w'] = None
            merged['expected_power_w'] = None
            merged['raw_pr'] = None

        if 'TIMESTAMP' in merged.columns:
            merged = merged.drop(columns=['TIMESTAMP'])
        per_group_dfs[group_id] = merged.reset_index()

    # Fixed MPPT (mppt 1) -- no tracker group
    fixed_df = meteo_data.copy()
    fixed_df['TIMESTAMP'] = pd.to_datetime(fixed_df['TIMESTAMP'])
    fixed_df = fixed_df.set_index('TIMESTAMP').sort_index()

    solar_pos = location.Location(-27.597, -48.549, tz='America/Sao_Paulo', altitude=3, name='Florianopolis').get_solarposition(fixed_df.index)
    dni_dhi = irradiance.erbs(ghi=fixed_df['GHIA_SMP22_Comp_Avg'], zenith=solar_pos['zenith'], datetime_or_doy=fixed_df.index)
    fixed_df['dni'] = dni_dhi['dni']
    fixed_df['dhi'] = dni_dhi['dhi']
    fixed_df.loc[fixed_df['GHIA_SMP22_Comp_Avg'] < 20, 'GHIA_SMP22_Comp_Avg'] = 0
    fixed_df['dni'] = fixed_df['dni'].where(fixed_df['GHIA_SMP22_Comp_Avg'] > 0, 0)
    fixed_df['dhi'] = fixed_df['dhi'].where(fixed_df['GHIA_SMP22_Comp_Avg'] > 0, 0)

    poa = irradiance.get_total_irradiance(
        surface_tilt=fixed_tilt,
        surface_azimuth=fixed_azimuth,
        solar_zenith=solar_pos['zenith'],
        solar_azimuth=solar_pos['azimuth'],
        dni=fixed_df['dni'],
        ghi=fixed_df['GHIA_SMP22_Comp_Avg'],
        dhi=fixed_df['dhi'],
        dni_extra=irradiance.get_extra_radiation(fixed_df.index),
        model='perez',
    )
    fixed_df['gpoa'] = poa['poa_global']

    # MPPT1 strings
    strings = mppt_string_map.get(fixed_mppt_id, [])
    cols = [f"dcw_{i}" for i in strings if f"dcw_{i}" in inverter_data.columns]
    if cols:
        fixed_df = fixed_df.join(inverter_data[cols], how='left')
        fixed_df['power_w'] = fixed_df[cols].sum(axis=1)
    else:
        fixed_df['power_w'] = None

    module_count = inverter_dict['modules_per_mppt'][fixed_mppt_id - 1]
    module_area_m2 = inverter_dict['module_area'][fixed_mppt_id - 1]
    area_m2 = module_count * module_area_m2
    fixed_df['area_m2'] = area_m2

    fixed_df['mppt_nominal_w'] = inverter_dict['mppts_power'][fixed_mppt_id - 1]
    irradiance_pu = fixed_df['gpoa'] / 1000
    fixed_df['expected_power_w'] = irradiance_pu * fixed_df['mppt_nominal_w']

    mask = (fixed_df['GHIA_SMP22_Comp_Avg'] >= 50) & (fixed_df['expected_power_w'] > 0)
    fixed_df['raw_pr'] = None
    fixed_df.loc[mask, 'raw_pr'] = fixed_df.loc[mask, 'power_w'] / fixed_df.loc[mask, 'expected_power_w']

    fixed_df['mppt_id'] = fixed_mppt_id
    fixed_df['group_id'] = None
    fixed_df['btactive_a1'] = False
    fixed_df['position_a1_degree'] = fixed_tilt
    fixed_df['targetangle_a1_degree'] = fixed_tilt

    # Build final per-MPPT dict of DataFrames
    mppt_dfs = {fixed_mppt_id: fixed_df.reset_index()}
    for group_id, df in per_group_dfs.items():
        mppt_id = group_to_mppt.get(int(group_id))
        if mppt_id not in mppt_dfs:
            mppt_dfs[mppt_id] = df.copy()
        else:
            mppt_dfs[mppt_id] = pd.concat([mppt_dfs[mppt_id], df], ignore_index=True, sort=False)

    # Rename columns for standardized headers
    for mppt_id, df in mppt_dfs.items():
        df.rename(columns={
            'GHIA_SMP22_Comp_Avg': 'ghi'
        }, inplace=True)
        df['target_position'] = df.get('targetangle_a1_degree', df['position_a1_degree'])
        # Remove individual dcw columns, keep only power_w
        dcw_cols = [col for col in df.columns if col.startswith('dcw_')]
        df.drop(columns=dcw_cols, inplace=True, errors='ignore')

    logging.info(f"Processed {len(mppt_dfs)} MPPTs for day {start_date}")
    return mppt_dfs


def main():
    parser = argparse.ArgumentParser(description='Calculate raw PR for AGRIPV MPPTs')
    parser.add_argument('start_date', help='Start date in YYYY-MM-DD format')
    parser.add_argument('end_date', help='End date in YYYY-MM-DD format')
    parser.add_argument('--csv', action='store_true', help='Export to CSV instead of uploading to DB')
    parser.add_argument('-d', '--debug', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO', help='Set logging level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.debug), format='%(asctime)s - %(levelname)s - %(message)s')

    start_dt = datetime.strptime(args.start_date, '%Y-%m-%d')
    end_dt = datetime.strptime(args.end_date, '%Y-%m-%d')

    date_range = pd.date_range(start=start_dt, end=end_dt, freq='D')
    for single_date in date_range:
        day_start = single_date.strftime('%Y-%m-%d')
        day_end = (single_date + timedelta(days=1)).strftime('%Y-%m-%d')
        day_mppt_dfs = process_day(day_start, day_end)
        if args.csv:
            out_dir = Path(__file__).resolve().parent
            for mppt_id, df in day_mppt_dfs.items():
                out_path = out_dir / f'mppt_{mppt_id}_{day_start}.csv'
                df.to_csv(out_path, index=False)
                logging.info(f"Exported MPPT {mppt_id} for {day_start} to {out_path}")
        else:
            inspector = inspect(engine_clean)
            for mppt_id, df in day_mppt_dfs.items():
                table_name = f'AGRIPV_raw_pr_mppt_{mppt_id:02d}'
                if inspector.has_table(table_name):
                    # Get existing columns
                    existing_columns = [col['name'] for col in inspector.get_columns(table_name)]
                    df_columns = df.columns.tolist()
                    missing_columns = [col for col in df_columns if col not in existing_columns]
                    if missing_columns:
                        # Add missing columns
                        with engine_clean.connect() as conn:
                            for col in missing_columns:
                                dtype = df[col].dtype
                                if pd.api.types.is_datetime64_any_dtype(dtype):
                                    sql_type = 'TIMESTAMP'
                                elif pd.api.types.is_bool_dtype(dtype):
                                    sql_type = 'BOOLEAN'
                                elif pd.api.types.is_numeric_dtype(dtype):
                                    sql_type = 'FLOAT'
                                else:
                                    sql_type = 'TEXT'
                                alter_sql = f'ALTER TABLE "{table_name}" ADD COLUMN "{col}" {sql_type}'
                                conn.execute(text(alter_sql))
                                conn.commit()
                # Now append (will create if not exists)
                df.to_sql(table_name, engine_clean, if_exists='append', index=False)
                logging.info(f"Uploaded MPPT {mppt_id} for {day_start} to table {table_name}")


if __name__ == "__main__":
    main()
