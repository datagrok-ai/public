# Row Source tests

**Filter values:**
- demog: `${AGE} > 44`
- spgi-100: `${Stereo Category} in ["R_ONE", "S_UNKN"]`

**Column setup reference:**

| Viewer | demog | spgi-100 |
|---|---|---|
| Scatter Plot | X = AGE, Y = HEIGHT, Color = RACE | X = Chemical Space X, Y = Chemical Space Y, Color = Stereo Category |
| Line Chart | X = AGE, Y = HEIGHT | X = Chemical Space X, Y = TPSA |
| Histogram | Column = AGE | Column = TPSA |
| Bar Chart | Value = AGE, Split = RACE | Value = TPSA, Split = Stereo Category |
| Pie Chart | Category = RACE | Category = Stereo Category |
| Box Plot | Category = RACE, Value = AGE | Category = Stereo Category, Value = TPSA |
| PC Plot | Columns: AGE, HEIGHT, WEIGHT | Columns: Chemical Space X, Chemical Space Y, TPSA |

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Open spgi-100

---

## Scatter Plot

1. Add Scatter Plot to demog (X = AGE, Y = HEIGHT, Color = RACE). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to M) — viewer should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — viewer should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — viewer should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — viewer should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — viewer should be empty
   10. Add a Pie Chart (Category = RACE) and hover over its slices — scatter plot should respond to hovering within the Filter setting; close the Pie Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update X to Chemical Space X, Y to Chemical Space Y, Color to Stereo Category; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Pie Chart (Category = Stereo Category) instead
6. Close the Scatter Plot

## Line Chart

1. Add Line Chart to demog (X = AGE, Y = HEIGHT). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to F) — viewer should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — viewer should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — viewer should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — viewer should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — viewer should be empty
   10. Add a Pie Chart (Category = RACE) and hover over its slices — line chart should respond to hovering within the Filter setting; close the Pie Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update X to Chemical Space X, Y to TPSA; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Pie Chart (Category = Stereo Category) instead
6. Close the Line Chart

## Histogram

1. Add Histogram to demog (Column = AGE). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to M) — histogram should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — histogram should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — histogram should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — histogram should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — histogram should be empty
   10. Add a Pie Chart (Category = RACE) and hover over its slices — histogram should respond to hovering within the Filter setting; close the Pie Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update Column to TPSA; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Pie Chart (Category = Stereo Category) instead
6. Close the Histogram

## Bar Chart

1. Add Bar Chart to demog (Value = AGE, Split = RACE). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to M) — bar chart should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — bar chart should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — bar chart should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — bar chart should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — bar chart should be empty
   10. Add a Pie Chart (Category = RACE) and hover over its slices — bar chart should respond to hovering within the Filter setting; close the Pie Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update Value to TPSA, Split to Stereo Category; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Scatter Plot (X = Chemical Space X, Y = Chemical Space Y) and lasso-select a group instead
6. Close the Bar Chart

## Pie Chart

1. Add Pie Chart to demog (Category = RACE). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to M) — pie chart should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — pie chart should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — pie chart should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — pie chart should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — pie chart should be empty
   10. Add a Bar Chart (Split = RACE) and hover over its bars — pie chart should respond to hovering within the Filter setting; close the Bar Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update Category to Stereo Category; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Bar Chart (Split = Stereo Category) instead
6. Close the Pie Chart

## Box Plot

1. Add Box Plot to demog (Category = RACE, Value = AGE). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to F) — box plot should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — box plot should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — box plot should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — box plot should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — box plot should be empty
   10. Add a Pie Chart (Category = RACE) and hover over its slices — box plot should respond to hovering within the Filter setting; close the Pie Chart
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update Category to Stereo Category, Value to TPSA; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Pie Chart (Category = Stereo Category) instead
6. Close the Box Plot

## PC Plot

1. Add PC Plot to demog (Columns: AGE, HEIGHT, WEIGHT). By default Row Source is set to **Filtered**
2. Set Filter to `${AGE} > 44` on the Context Panel
3. Check Row Sources:
   1. Add a filter on the Filter Panel (e.g., restrict SEX to M) — PC Plot should respond to both the Filter Panel and the Filter setting simultaneously
   2. Set Row Source to **All** — PC Plot should respond to the Filter setting only; remove Filter Panel filter
   3. Set Row Source to **Selected** — PC Plot should be empty
   4. Reset Filter Panel filter; select some rows — only selected rows where AGE > 44 should be displayed
   5. Set Row Source to **SelectedOrCurrent** — the same selected rows should be displayed
   6. Reset selection and make any row current — only the current row should be displayed (considering the Filter setting)
   7. Set Row Source to **FilteredSelected** — PC Plot should be empty
   8. Select rows where AGE is 42–47: only rows where AGE is 45, 46, or 47 should be displayed
   9. Set Row Source to **MouseOverGroup** — PC Plot should be empty
   10. Add a Histogram (Column = AGE) and hover over its bins — PC Plot should respond to hovering within the Filter setting; close the Histogram
   11. Set Row Source to **CurrentRow** — only the current row where AGE > 44 should be displayed
   12. Set Row Source to **MouseOverRow** — only the hovered row where AGE > 44 should be displayed
4. Set Table to spgi-100 (Table (2)) in Context Panel > Data; update Columns to Chemical Space X, Chemical Space Y, TPSA; set Row Source to **Filtered**; set Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Repeat step 3 using the spgi-100 data; for the Filter Panel filter add a TPSA range filter; for FilteredSelected: select rows spanning multiple Stereo Category values — only R_ONE/S_UNKN rows from the selection should be displayed; for MouseOverGroup: add a Histogram (Column = TPSA) instead
6. Close the PC Plot

---

7. Close demog
8. Close spgi-100

---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/spgi-100.csv"]
}
