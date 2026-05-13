### Principal component analysis (PCA)

1. Open the DataFrame:
* Open the cars.csv file from the Demo files.
* Alternatively, click on the 'star' icon in TestTrack (tooltip: "Open test data") to load the cars.csv dataset.
2. Run Principal Component Analysis (PCA): Go to Top Menu > ML > Analyze > PCA….
3. Select Features: In the Features input box, select all available columns in the dataset. Set Components to 3.
4. Execute PCA by clicking OK.
* Expected Results: Verify that three new columns (PC1, PC2, PC3) are added to the dataset.

5. Repeat PCA with Center and Scale Options:

Go to Top Menu > ML > Analyze > PCA…. In the Features input box, select all available columns in the dataset. Set Components to 3. Check the Center and Scale checkboxes. Click OK.
* Expected Results: Verify that three additional columns (PC1(2), PC2(2), PC3(2)) are added to the dataset.

---
{
"order": 1,
"datasets": ["System:DemoFiles/cars.csv"]
}