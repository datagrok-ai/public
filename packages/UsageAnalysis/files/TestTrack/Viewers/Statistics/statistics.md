---
feature: statsviewer
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.stats-viewer, entities.viewer.action.open-as-table, entities.viewer.action.full-screen]
realized_as:
  - statistics-spec.ts
related_bugs: []
---

# Statistics viewer (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Statistics** from **Toolbox > Viewers**

## Add, close and re-add

1. Click the **Statistics** icon in **Toolbox > Viewers** — the viewer opens with
   statistics drawn
2. Close it from the title bar — the viewer is gone
3. Click the icon again — it opens and draws again

## Adding a statistic

1. Right-click the viewer and open the **Statistics** submenu
2. `sum` is unticked
3. Click `sum` — the viewer redraws with the extra column

   Note: the submenu opens only once per page. Reopening it (to remove `sum`
   again, or to reach **Histograms**) does not work under automation — those
   steps are manual, below.

## Row source

1. On **Context Panel > Data**, set **Row Source** to *Filtered*
2. Filter the table to `SEX = M` — the statistics are recomputed and redraw
3. Reset the filter, set **Row Source** to *Selected*
4. Select every second row — the statistics redraw
5. Set **Row Source** back to *All*

## Columns

1. **Context Panel > Data > Columns** reports as many columns as the table has

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Manual scenarios (not automated)

Everything below was in the original checklist and is **not** covered by
`statistics-spec.ts`. Kept verbatim so no scenario is lost.

### Default statistics display

> Manual

1. Verify the viewer shows stat columns: values, nulls, unique, min, max, avg, med, stdev
2. Verify numerical columns AGE, HEIGHT, WEIGHT have non-empty values in all stat columns
3. Verify the name column lists the demog dataset columns

### Statistics for categorical columns

> Manual

1. Verify numerical stat columns (avg, min, max, stdev) are empty for categorical column rows (SEX, RACE, DIS_POP)
2. Verify count-type stat columns (values, nulls, unique) are populated for categorical column rows

### Add and remove statistics columns

> Manual

1. Right-click the Statistics viewer
2. Hover over Statistics in the context menu — submenu opens listing available stat types
3. Click sum in the submenu — a sum column appears in the viewer
4. Right-click the viewer and hover over Statistics in the context menu
5. Click sum again — the sum column disappears

### Histogram columns

> Manual

1. Right-click the Statistics viewer
2. Hover over Histograms in the context menu — submenu lists categorical columns with fewer than 10 categories (SEX, RACE, DIS_POP)
3. Click SEX in the submenu — a histogram sparkline column appears in the viewer
4. Right-click the viewer and hover over Histograms in the context menu
5. Click SEX again — the histogram column disappears

### Open as table

> Manual

1. Right-click the Statistics viewer
2. Click Open as table — a new table tab opens with stat types as column headers and data columns as rows

### Full-screen shortcut

> Manual

1. Click the expand icon on the Statistics viewer panel — viewer expands to full screen
2. Press Alt+F — viewer returns to normal size

---
{
  "order": 100,
  "datasets": ["System:DemoFiles/demog.csv"]
}
