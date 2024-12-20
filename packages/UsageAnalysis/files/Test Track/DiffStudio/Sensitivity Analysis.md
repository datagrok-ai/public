### Sensitivity Analysis in Diff Studio (Bioreactor)

To validate the functionality of the Sensitivity Analysis tool in [Diff Studio](https://docs.google.com/document/d/1vXeUMv6S8L94KvMsOHgodOiJ2eHzsuHRzoFGnONm0jg/edit?tab=t.0), ensuring that the analysis runs correctly and the results are visualized properly.

Test Steps:

1. Open Diff Studio and Load Bioreactor Model:

- Go to Apps and run Diff Studio > Click on the Open model icon on the ribbon > Select Examples > Bioreactor to load the model.
- Turn off the *Edit* toggle on the ribbon; equations editor closes, form with inputs opens

2. Open Sensitivity Analysis:

Click on the *Sensitivity* icon on the ribbon. The Sensitivity Analysis view should open.

3. Select Parameters for Analysis:

In the Sensitivity Analysis view, use the switchers to select the following parameters: modify *Process mode*; check that *FFox & KKox* (and some other) inputs are modified.

4. Run Sensitivity Analysis. Click the Run icon on the ribbon to start the analysis.

Observe the Results: After running the analysis, four viewers should open and visualize the results.

**Expected Results:**

* The Sensitivity Analysis view should open without errors.
* The parameters (FFox, FKox, FFred) should be selectable without any issues.
* Four viewers should open, each displaying the sensitivity analysis results. The visualizations should be clear and accurate, with no errors.

---
{
  "order": 3
}