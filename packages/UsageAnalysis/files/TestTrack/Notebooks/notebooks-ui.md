---
feature: notebooks
realizes_atlas: [notebooks.cp.browse-and-open-html]
realizes: [notebooks.view.notebook]
priority: p0
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  The JupyterLab editor runs inside a sandboxed iframe that automated browser
  tests cannot reach, so the edit-and-run steps must be run by hand.
related_bugs:
  - GROK-6630
---




## Scenarios

### Create Notebook

1.Open "demog" table.
2. Right-click the table view title and select Open in Notebook from the context menu
3. On the top toolbar, click the HTML button.
4. On the top toolbar, click the Download button and download the created notebook in all suggested formats.


### Edit and run a notebook


1. Go to **Browse > Platform > Notebooks**.
2. Find the notebook named "demog" and double-click it.
   - Expected: the notebook opens in read-only HTML view.
3. In the menu ribbon, click **Edit**.
   - Expected: the notebook opens in JupyterLab edit mode inside the editor iframe.
4. Add the following code to the notebook body:

   ```python
   race = demog['RACE']
   race.describe()
   ```

   - Expected: the code is visible in the active notebook cell.
5. In the menu ribbon, click **Play**.
   - Expected: the cell executes and the output (a summary of the "RACE" column) is displayed
     below the cell with no kernel errors.



## Cleanup

- Delete the demog notebook.

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv"]
}
