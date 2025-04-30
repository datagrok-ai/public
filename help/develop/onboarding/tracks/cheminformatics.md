---
title: "Cheminformatics"
---

These programming exercises are designed to introduce developers to the Datagrok platform cheminformatics capabilities.
The exercises are based on your knowledge obtained in [exercises](../exercises.md).

## Table of contents

* [Search for most common structures](#exercise-1-search-for-most-common-structures)
* [Model to predict activity](#exercise-3-train-model-to-predict-activity)

## Basic exercises in cheminformatics

### Exercise 1: Search for most common structures

_You will learn:_ How to employ functions from external packages in your own package.

_Prerequisites:_ ["Cheminformatics"](../../../datagrok/solutions/domains/chem/chem.md).

_Statement of the problem_. Write a function that reads a file containing SMILES, determines the associated maximal
common substructure (MCS), and computes the mutual similarity scores for molecules and the MCS.

_Input data_. Files > App Data > Chem > sars\_small.csv

**Solution**, step-by-step.

1. Let's call our function `findSimilarToMCS`, we place its definition in `./src/package.ts` inside our package. This
   function takes a single input — a dataframe `df`. For the sake of simplicity, we suppose that the column with SMILES
   is `df.col('smiles')`:

    ```typescript
    //name: findSimilarToMCS
    //input: dataframe df
    export async function findSimilarToMCS(df: DG.DataFrame) : Promise<void> {
      ... // your code goes here
    }
    ```

2. Employ the asynchronous function `FindMCS` from `Chem` package. Since we're calling a function from an external
   package, we should use `grok.functions.call`:

    ```typescript
    const mcsValue = await grok.functions.call('Chem:FindMCS', {'smiles': 'smiles', 'df': dataframe, 'returnSmarts':
    false});
    ```

3. Having obtained the string `mcsValue`, create a new column in `df`, whose cells are filled with this value:
     * Create an `Array` of the appropriate length, filled with `mcsValue`.
     * Feed this array to the constructor `DG.Column.fromList()` to get the desired `mcsCol` object.
     * Assign semantic type `Molecule` to the newly created column, with the help of `col.semType(...)`. Similarly,
       associate `Molecule` cell renderer with the help of `col.setTag(...)` method.

4. To compute similarity scores, we can call the `getSimilarities()` function of `Chem` package, which takes as its
   parameters the initial SMILES column and `mcsValue`. The function can be invoked as described in step 2.

5. The output of step 4 is a new dataframe `scoresDf`, its 0-th column contains the scores values. This
   column, `scoresCol`, can be reached with by means of `byIndex()` method of `scoresDf.columns` object.

6. Finally, insert the columns `mcsCol` and `scoresCol` into the dataframe, next to the position of the initial SMILES
   column. `df.columns.insert()` method can help with this, if we cleverly specify the index/position at which the
   insertion should take place.

### Exercise 3: Train Model to Predict Activity

_You will learn:_ How to train a model inside a package and use it to predict the activity of molecules

_Prerequisites:_  ["Molecular fingerprints"](../../../datagrok/solutions/domains/chem/fingerprints.md),
["Cheminformatics"](../../../datagrok/solutions/domains/chem/chem.md).

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

3. Using [grok.chem.descriptors](https://datagrok.ai/js-api/dg/namespaces/chem/functions/descriptors)  create fingerprint of all
   molecules.
4. Use grok.ml.trainModel your model (using fingerprint) to predict activity of molecule. You can use
   dataset [example](https://public.datagrok.ai/f/Demo.TestJobs.Files.DemoFiles/chem/activity_cliffs.csv)
5. Using [grok.ml.applyModel](https://datagrok.ai/js-api/dg/namespaces/ml/functions/applyModel) apply on the test and train
   datasets. Check the accuracy of the model.
6. Using grok.shell.addTableView(datasetName) output test dataset
