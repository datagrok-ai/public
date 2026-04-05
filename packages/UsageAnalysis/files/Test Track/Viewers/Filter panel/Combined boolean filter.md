1. Open **demog**
2. Open the **Filter Panel**
3. Add a new column: **SEX_bool**: `case ${SEX} when "F" then true else false end`
5. Verify that each boolean column (CONTROL, SEX_bool) appears as a row in the Combined Boolean filter with the label "Flags OR", and count indicators are shown (CONTROL = 39, SEX_bool = 3243).
7. Apply the combined filter: for CONTROL check the first checkbox (true) and uncheck the second (false); for SEX_bool uncheck the first (true) and check the second (false). In OR mode this filters to CONTROL=true OR SEX_bool=false — verify filtered row count = 2632.
8. Apply other filters: in RACE categorical filter check only `Asian` (other RACE values show 0). In AGE histogram filter set range to 50–89.
9. Save the layout via Toolbox → Layouts → SAVE.
10. Close the Filter Panel
11. Apply the saved layout by clicking the layout thumbnail — verify the filter state and row count match the state before close.
12. Remove all filters via Hamburger menu → Remove All.
13. Close the **Filter Panel**.
14. Open the **Filter Panel** — the Combined Boolean filter should be added automatically.
15. Apply the saved layout by clicking the layout thumbnail — verify the filter state and row count match the state before close.

---
{
"order": 8,
"datasets": ["System:DemoFiles/demog.csv"]
}
