# PMX Studio demo data

Studies to try PMX Studio with. The quickest way: on the PMX Studio Overview, each stage's card has **Try with demo
data** at its foot (Plan, Monitor and Analyze open with the study below), and while you have no studies yet, **Or
compare two demo studies** saves the CT26 and MC38 studies for you and opens them in Compare.

The files themselves are on the server in **Browse > Files > App Data > Pmxstudio > demo**. Open one from the Overview
(**Open a study**) or a stage's **Open…**, or drag it onto the app.

All are synthetic (made by seeded models: not real animals, not real PK), except the published MC38 study in
`../example/` (Yu X. et al., *Nat Commun* 15 (2024), doi:10.1038/s41467-024-55006-2, CC BY 4.0).

## What to open

| To try | Open |
|---|---|
| Plan a study | `ct26-trametinib-antipd1.pmxplan` |
| A model the atlas has not seen | `PDX-SYN-01 pilot plan.pmxplan`, or in Plan, **+ Add your own model** with `synthetic-history-PDX-SYN-01.csv` |
| Monitor a running study | `ct26-running-day9.pmxstudy`, or `ct26-trametinib-antipd1-running-day9.csv` with the plan |
| Analyze a finished study | `ct26-trametinib-antipd1-finished.csv` (with the plan), or `../example/mc38-ld4172-anti-pd1.csv` (a combination) |
| Exposure against patients (needs the NCA package) | a finished study, then its `*-pk-mouse.csv` in Analyze's exposure |
| Compare studies | save `ct26-trametinib-antipd1-finished.csv` and `mc38-trametinib-finished.csv` with their PK, tick both in the Overview |
| Other table layouts | `mc38-example-wide-by-day.csv` (a column per day), `mc38-example-calipers-dates.csv` (calipers and dates), `xlsx-groups-01-and-1.xlsx` (Excel, group names that look alike) |

## Files

- `ct26-trametinib-antipd1-finished.csv`: CT26, 6 groups × 10 mice (vehicle; trametinib 0.3, 1, 3 mg/kg p.o. QD;
  anti-PD-1 10 mg/kg i.p. BIW; trametinib 1 + anti-PD-1), days 0–21, body and tumor weight, removal reasons. The
  3 mg/kg group crosses the welfare limits on purpose.
- `ct26-trametinib-antipd1-running-day9.csv`: the same design measured to day 9, with planted cases: an animal over
  the 1500 mm³ limit, one growing fast, two losing weight, a caliper typo ×10, a missed measurement.
- `ct26-trametinib-antipd1.pmxplan`, `ct26-running-day9.pmxstudy`, `PDX-SYN-01 pilot plan.pmxplan`: saved by PMX Studio
  from these files.
- `ct26-trametinib-pk-mouse.csv`, `mc38-trametinib-pk-mouse.csv`: mouse plasma PK after the first dose (day 0), 4 animals per
  trametinib group (columns SUBJECT, GROUP, COMPOUND, ROUTE, DOSE_MG_KG, TIME_H, CONC_NG_ML).
- `mc38-trametinib-finished.csv`: a second, MC38 trametinib study (vehicle, 1 and 3 mg/kg), days 0–18.
- `synthetic-history-PDX-SYN-01.csv`: an implantation pilot of a PDX outside the atlas, 2 groups × 6 mice, days 0–28.
- `mc38-example-*.csv`: the published MC38 study in other layouts. `xlsx-groups-01-and-1.xlsx`: an Excel table whose
  group names "01" and "1" must stay apart.

The NCA package's own demo studies are in **App Data > Nca**.
