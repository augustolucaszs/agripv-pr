# AGRIPV — Performance Ratio Calculation

Daily Performance Ratio (PR) pipeline for the agrivoltaic (AGRIPV) system located at the
Federal University of Santa Catarina (UFSC), Florianópolis, Brazil (−27.597°S, −48.549°W).

---

## System Description

The AGRIPV system combines solar PV panels mounted on single-axis trackers over an agricultural
area. It consists of 7 MPPT channels fed by one SMA inverter:

| MPPT | Nominal Power | Type            | DC String | Tracker Groups |
|------|--------------|-----------------|-----------|----------------|
| 1    | 10.695 kW    | Fixed (N-facing, 15°) | dcw\_1 | — |
| 2    | 9.265 kW     | Single-axis tracker   | dcw\_2 | Group 2 |
| 3    | 8.720 kW     | Single-axis tracker   | dcw\_3 | Group 3 |
| 4    | 8.720 kW     | Single-axis tracker   | dcw\_4 | Group 4 |
| 5    | 9.265 kW     | Single-axis tracker   | dcw\_5 | Group 5 |
| 6    | 12.250 kW    | Single-axis tracker   | dcw\_6 | Groups 6 & 7 |
| 7    | 12.250 kW    | Single-axis tracker   | dcw\_7 | Groups 8 & 9 |

Trackers are north–south horizontal single-axis (HSAT). The tracker angle convention is:
negative = east tilt (morning), positive = west tilt (afternoon), 0° = horizontal.

MPPTs 6 and 7 each span two physical tracker rows. Their tracker positions are averaged to
obtain a representative tilt angle for irradiance calculation.

---

## Data Sources (PostgreSQL — `fotovoltaica` database)

| Table | Description | Key columns |
|-------|-------------|-------------|
| `agripv_inverter_sma` | 1-min DC power per string | `TIMESTAMP`, `dcw_1`–`dcw_12` (W) |
| `SapAlbedo_1m` | 1-min meteorological data | `TIMESTAMP`, `GHIA_SMP22_Comp_Avg` (GHI, W/m²) |
| `trackers_suntrack_tcu` | 1-min tracker telemetry | `TIMESTAMP`, `group_id`, `position_a1_degree`, `targetangle_a1_degree`, `btactive_a1` |

---

## Methodology

### Plane-of-Array (POA) Irradiance

Since no direct POA pyranometer is mounted on each tracker, the effective irradiance on each
MPPT's panel surface is computed from the measured Global Horizontal Irradiance (GHI) using
**pvlib**:

1. **Erbs decomposition** — splits GHI into Direct Normal Irradiance (DNI) and Diffuse
   Horizontal Irradiance (DHI).

2. **Perez sky diffuse model** — accounts for circumsolar and horizon brightening to compute
   the total POA irradiance (`poa_global`) at each tracker tilt angle:

   ```
   GPOA = GHI_beam × cos(AOI) + DHI_Perez + GHI × albedo × (1 − cos(tilt)) / 2
   ```

3. **Surface orientation per minute:**
   - *Tracked MPPTs (2–7):* `surface_tilt = |position_a1_degree|`, `surface_azimuth = solar_azimuth`
     (N–S axis → panel azimuth tracks the sun east–west).
   - *Fixed MPPT 1:* `surface_tilt = 15°`, `surface_azimuth = 0°` (north-facing).

4. Low-irradiance filter: GHI values below 20 W/m² are set to zero before decomposition to
   avoid numerical artefacts from Erbs at near-zero irradiance.

### Performance Ratio (IEC 61724-1)

The PR is computed as the ratio of daily energy sums, **not** a mean of per-minute ratios:

```
         Σ [ P_actual(t) · Δt ]
PR_day = ─────────────────────────────────────────
         Σ [ GPOA(t) · Δt ] / 1000 × P_nom
```

Where:
- `P_actual(t)` — measured DC string power at minute `t` (W)
- `GPOA(t)` — computed plane-of-array irradiance at minute `t` (W/m²)
- `Δt` — 1 minute = 1/60 h (cancels in numerator and denominator)
- `P_nom` — MPPT nominal DC power (W)
- Division by 1000 converts irradiance from W/m² to kW/m² for unit consistency

Only minutes where **GHI ≥ 50 W/m²** are included in the sums. This excludes dawn/dusk
periods where inverter measurement lag and diffuse scatter cause unreliable per-minute ratios.

The formula simplifies to:

```
         Σ P_actual(t)
PR_day = ──────────────────────────
         Σ GPOA(t) / 1000 × P_nom
```

A PR of 1.0 (100%) means the system produced exactly as much energy as a reference system
with the same nominal power operating at STC efficiency under the actual POA irradiation.
Values above 1.0 are possible due to spectral effects, cooler temperatures (modules perform
better below 25 °C), or minor irradiance model errors on partly cloudy days.

---

## Script — `AGRIPV_raw_pr.py`

### Requirements

```
pvlib
pandas
numpy
sqlalchemy
psycopg2
python-dotenv
```

### Configuration

Copy `.env.example` to `.env` and fill in database credentials:

```
DB_SOURCE_URL=postgresql://user:password@host/fotovoltaica
DB_CLEAN_URL=postgresql://user:password@host/fotovoltaica_clean
```

The `.env` file is listed in `.gitignore` and must never be committed.

### Usage

```bash
# Process a date range and push to the clean DB
python3 AGRIPV_raw_pr.py 2025-09-08 2026-06-15

# Test a single day — write CSV files instead of pushing to DB
python3 AGRIPV_raw_pr.py 2025-09-09 2025-09-09 --csv

# Push only the daily PR summary table (skip per-minute tables)
python3 AGRIPV_raw_pr.py 2025-09-08 2026-06-15 --daily-only

# Verbose debug output
python3 AGRIPV_raw_pr.py 2025-09-09 2025-09-09 --csv -d DEBUG
```

### Processing Flow (per day)

```
Source DB ──► Fetch inverter, meteo, tracker data
               │
               ▼
         Prep & align timestamps
               │
               ├─► Fixed MPPT 1
               │     tilt=15°, azimuth=0° (fixed)
               │     GHI → Erbs → Perez → GPOA
               │
               └─► Tracked MPPTs 2–7
                     Average tracker position across groups (for split MPPTs 6 & 7)
                     GHI → Erbs → Perez → GPOA using |tracker angle|
                     Join DC power from corresponding dcw_N string
               │
               ▼
         Per-minute DataFrame per MPPT
         [TIMESTAMP, mppt_id, power_w, gpoa, ghi, position_deg, target_deg, tracker_active]
               │
               ├─► Upload to AGRIPV_raw_pr_mppt_01..07  (per-minute, clean DB)
               │
               ▼
         Compute daily PR (IEC 61724) per MPPT
         Filter: GHI ≥ 50 W/m²
         PR = Σ(power_w) / (Σ(gpoa) / 1000 × nominal_kw)
               │
               └─► Upload to AGRIPV_daily_pr  (one row per MPPT per day, clean DB)
```

---

## Output Tables (`fotovoltaica_clean` database)

### `AGRIPV_raw_pr_mppt_01` through `AGRIPV_raw_pr_mppt_07`

One row per minute per MPPT.

| Column | Type | Description |
|--------|------|-------------|
| `TIMESTAMP` | timestamp | 1-minute measurement time |
| `mppt_id` | int | MPPT index (1–7) |
| `power_w` | float | DC string power (W) |
| `gpoa` | float | Computed plane-of-array irradiance (W/m²) |
| `ghi` | float | Global horizontal irradiance — clipped (W/m²) |
| `position_deg` | float | Tracker tilt angle (°); negative = east, positive = west |
| `target_deg` | float | Tracker target angle (°) |
| `tracker_active` | bool | Whether tracker control was active |

### `AGRIPV_daily_pr`

One row per MPPT per day.

| Column | Type | Description |
|--------|------|-------------|
| `date` | text | Date (YYYY-MM-DD) |
| `mppt_id` | int | MPPT index (1–7) |
| `pr` | float | IEC 61724 daily Performance Ratio |
| `energy_kwh` | float | DC energy produced during GHI ≥ 50 W/m² period (kWh) |
| `poa_kwh_m2` | float | POA irradiation during same period (kWh/m²) |
| `nominal_kw` | float | MPPT nominal power (kW) |
| `n_minutes` | int | Number of 1-min data points used in PR calculation |
| `tracking` | bool | Whether MPPT uses a tracker |

---

## Grafana Queries

### Daily PR — all MPPTs in one panel

```sql
SELECT
  to_timestamp("date", 'YYYY-MM-DD')         AS time,
  pr                                           AS value,
  'MPPT ' || lpad(mppt_id::text, 2, '0')     AS metric
FROM "AGRIPV_daily_pr"
WHERE "date"::date BETWEEN $__timeFrom()::date
                       AND $__timeTo()::date
  AND pr IS NOT NULL
  AND pr < 1.5
ORDER BY time, mppt_id
```

### Instantaneous DC power ratio — all MPPTs (per-minute)

```sql
SELECT "TIMESTAMP" AS time, power_w/(gpoa/1000.0*10.695) AS value, 'MPPT 01' AS metric FROM "AGRIPV_raw_pr_mppt_01" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0 AND power_w/(gpoa/1000.0*10.695)<1.5
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*9.265),  'MPPT 02' FROM "AGRIPV_raw_pr_mppt_02" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*8.720),  'MPPT 03' FROM "AGRIPV_raw_pr_mppt_03" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*8.720),  'MPPT 04' FROM "AGRIPV_raw_pr_mppt_04" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*9.265),  'MPPT 05' FROM "AGRIPV_raw_pr_mppt_05" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*12.250), 'MPPT 06' FROM "AGRIPV_raw_pr_mppt_06" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
UNION ALL
SELECT "TIMESTAMP", power_w/(gpoa/1000.0*12.250), 'MPPT 07' FROM "AGRIPV_raw_pr_mppt_07" WHERE $__timeFilter("TIMESTAMP") AND ghi>=50 AND gpoa>0
ORDER BY time, metric
```

Set **Format as: Time series** in the Grafana query editor for both queries.

---

## Known Limitations

- **MPPT2 underperformance:** MPPT2 shows PR ≈ 7–10% for most of the dataset due to a hardware
  fault. Exclude it from fleet-level averages: `WHERE mppt_id != 2`.

- **No direct POA sensor per tracker:** GPOA is modelled from GHI using Erbs + Perez. On
  partly cloudy days with rapid irradiance transients, the 1-minute averaging can produce
  instantaneous PR outliers. The daily aggregation (ratio-of-sums) is robust to this.

- **Bifacial rear-side irradiance:** Rear-side irradiance (APOA) is measured but not currently
  included in the PR reference irradiance. The current PR uses front-side GPOA only, so on
  good days the effective bifacial gain appears as PR > 1.

---

## References

- IEC 61724-1:2021 — *Photovoltaic system performance — Part 1: Monitoring*
- Anderson, K. et al. — pvlib python: a python package for modeling solar energy systems.
  *Journal of Open Source Software*, 2023.
- Erbs, D.G., Klein, S.A., Duffie, J.A. — Estimation of the diffuse radiation fraction for
  hourly, daily and monthly-average global radiation. *Solar Energy*, 1982.
- Perez, R. et al. — Modeling daylight availability and irradiance components from direct and
  global irradiance. *Solar Energy*, 1990.
