### Fitting in Diff Studio (Bioreactor)
To validate the functionality of the Fitting tool in Diff Studio, ensuring that the fitting process is executed correctly and the results are accurately visualized.

1. Open **Diff Studio** and Load **Bioreactor Model**:
* Go to Apps and run Diff Studio. Click on the Open model icon on the ribbon.
* Select Examples > Bioreactor to load the model.
2. Open **Fitting View**. Click on the Run fitting icon on the ribbon. The Fitting view should open without errors.
3. Select **Parameters for Fitting**:
* In the Fitting view, use the switchers to select the following parameters:
  * switch at
  * FFox from 0.15 to 1.0
  * FKox from 0 to 3
4. Scroll to the **Target Block** to locate the Target block in the Fitting view. 
5. Input **Bioreactor Data**: In the Bioreactor table input field, add the file from Files: App Data > Diff Studio > examples > bioreactor-experiment.csv
6. Run the Fitting Process. Click the **Run** icon on the ribbon to start the fitting process.

**Expected Results**: The parameters (switch at, FFox, FKox) should be selectable without any issues. The Bioreactor.csv file should be added successfully to the Bioreactor table input. The "RMSE by iterations" column should display a descending graph, reflecting the fitting results accurately. No errors should occur during the process.

---
{
  "order": 4
}