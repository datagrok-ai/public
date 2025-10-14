### Fitting in Diff Studio (Bioreactor)
To validate the functionality of the Fitting tool in Diff Studio, ensuring that the fitting process is executed correctly and the results are accurately visualized.

1. Open **Diff Studio** and Load **Bioreactor Model**:
* Go to Apps and run Diff Studio. Click on the Open model icon on the ribbon.
* Select Examples > Bioreactor to load the model.
2. Click on the Fit icon on the ribbon, the Fitting view opens without errors.
3. Modify *Process mode*; check that *FFox & KKox* (and some other) inputs are modified
4. Set *Process mode* to “Default”
* use the switchers to select the following parameters:
  * switch at
  * FFox from 0.15 to 1.0
  * FKox from 0 to 3
4. Scroll to the **Target Block** to locate the Target block in the Fitting view. 
5. Input **Bioreactor Data**: In the Bioreactor table input field, add the file from Files: App Data > Diff Studio > library > bioreactor-experiment.csv
6. Run the Fitting Process. Click the **Run** icon on the ribbon to start the fitting process.
(REMARK. Grid may contain another number of rows)


**Expected Results**: The parameters (switch at, FFox, FKox) should be selectable without any issues. The Bioreactor.csv file should be added successfully to the Bioreactor table input. The "RMSE by iterations" column should display a descending graph, reflecting the fitting results accurately. No errors should occur during the process.

---
{
  "order": 4
}