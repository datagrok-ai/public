<!-- TITLE: Cheminformatics -->
<!-- SUBTITLE: -->

# Cheminformatics

These programming exercises are designed to introduce developers to the Datagrok platform cheminformatics capabilities.
The exercises are based on your knowledge obtained in [exercies](./exercises.md).

## Table of contents

* [Search for most common structures](#exercise-1-search-for-most-common-structures)
* [Model to predict activity](#exercise-3-train-model-to-predict-activity)

# Basic exercises in cheminformatics

## Exercise 1: Search for most common structures

Statement of the problem.
Write a function that reads a file containing SMILES, determines the associated maximal common substructure (MCS).
It should also compute the mutual similarity scores for molecules and the MCS.
Input data. Files > App Data > Chem > sars\_small.csv
*Solution*, step-by-step.

1. Let's call our function findSimilarToMCS. This function takes a single input -- a dataframe df.
For the sake of simplicity, we suppose that the column with SMILES is df.col('smiles').
2. Employ the asynchronous function FindMCS from Chem package.
Since we're calling a function from an external package, we should use grok.functions.call:
const mcsValue = await grok.functions.call('Chem:FindMCS', {'smiles': 'smiles', 'df': dataframe, 'returnSmarts': false});
3. Having obtained the string mcsValue, create a new column in df, whose cells are filled with this value:
Create an Array of the apporpriate length, filled with mcsValue.
Feed this array to the constructor DG.Column.fromList() to get the desired mcsCol object.
Assign semantic type Molecule to the newly created column, with the help of col.semType(...).
  Similarly, associate Molecule cell renderer with the help of col.setTag(...) method.
4. To compute similarity scores, we can invoke the getSimilarities() function of Chem.
It takes as its parameters the initial SMILES column and mcsValue. The function can be called as described in step 2.
5. The output of step 4 is a new dataframe scoresDf, its 0-th column contains the scores values.
This column, scoresCol, can be reached with by means of byIndex() method of scoresDf.columns object.
6. Finally, insert mcsCol and scoresCol next to the position of the initial SMILES column.

## Exercise 3: Train Model to Predict Activity

*You will learn:* How to train a model inside a package and use it to predict the activity of molecules

*Prerequisites:*  ["Molecular fingerprints"](#https://datagrok.ai/help/domains/chem/fingerprints),
["Cheminformatics"](#https://datagrok.ai/help/domains/chem/cheminformatics).

1. Create a package with the name `<yourFirstName>-cheminformatics`
2. Add new function

   ```javascript
   // name: TrainAndPredict
   //input: dataframe train
   //input: dataframe test
   export function TrainAndPredict(train, test) {
      // your code here
   }
   ```

Here the training and test dataframes are our datasets for training and prediction, respectively.
3. Using [grok.chem.descriptors](#https://datagrok.ai/js-api/modules/grok.chem#descriptors)  create fingerprint of all molecules.
4. Use [grok.ml.trainModel](#unexist) your model (using fingerprint) to predict activity of molecule.
You can use dataset [example](#https://public.datagrok.ai/f/Demo.TestJobs.Files.DemoFiles/chem/activity_cliffs.csv)
5. Using [grok.ml.applyModel](#https://datagrok.ai/js-api/modules/grok.ml#applyModel) apply on the test and train datasets.
Check the accuracy of the model.
6. Using grok.shell.addTableView(datasetName) output test dataset
