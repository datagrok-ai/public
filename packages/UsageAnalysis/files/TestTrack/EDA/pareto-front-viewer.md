### Pareto front viewer test case

#### Empty & Non-Numeric columns handling
1. Open **cars-with-missing.csv** data from Demo files (or press Star icon). This table has an empty column - ''turbo'';
2. Add the Pareto front viewer. Open the viewer Properties panel.
3. In the Objectives section, check the **Minimize** and **Maximize** dropdown lists.
- Expected: Empty columns (e.g., turbo) must **NOT** be present in **Minimize** or **Maximize** lists.
- Expected: Only numeric columns should appear in these dropdowns — e.g. the model column (string) must **NOT** be included.
4. Select **all** columns in **Maximize**. 
- Expected: A warning must appear indicating that it is impossible to apply both minimization and maximization to the same feature simultaneously.

#### Category column & Label selection behavior

5. Open cars.csv from Demo files and add the Pareto Front viewer.
- Expected: The **model** column is automatically selected as **Label**, as it contains unique values.
6. Open the **demog** dataset (via Tools > …). Add the Pareto Front viewer.
- Expected: The Label column selector is empty by default.
- Expected: The viewer should auto-select a category column **only if it has unique values**. If no such column exists, Label remains empty. 

#### Properties & UI behavior validation
7. Review all viewer properties: Labels, Objectives, Axes.
- Expected: No UI errors, no missing options, no unexpected behavior.
- Expected: Changing objectives/labels/axes updates the viewer correctly.


---
{
"order": 5,
"datasets": ["System:DemoFiles/cars-with-missing.csv"]
}
