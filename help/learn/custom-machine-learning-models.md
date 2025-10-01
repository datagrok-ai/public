---
title: "Custom machine learning models"
---

Datagrok supports two ML engines out of the box: Caret and Chemprop. In addition, the platform allows users to build
their own customizable ML models. This capability provides full access and control to algorithms available in ML
libraries in any of the [supported languages](../compute/compute.md#functions-and-cross-language-support). The users may construct and
configure any chosen model in their custom data pipeline.

Custom models can be configured in several ways:  

* **Script-based functions** — straightforward Python, R, or Julia scripts embedded in Datagrok (see [example](#example) below).
* **Package functions** — functions provided as part of a Datagrok package. For example, see the [TensorFlow.js](https://github.com/datagrok-ai/public/tree/master/packages/TensorFlow.js) package.
* **Celery tasks inside Docker containers** — ideal for complex models. Benefits include:  
  * GPU support for heavy computations  
  * Isolated, reproducible environments  
  * Easier modularization and maintenance  

For examples and a detailed setup, see the [Samples package](https://github.com/datagrok-ai/public/tree/master/packages/Samples/dockerfiles/) and the [Celery Python functions documentation](https://datagrok.ai/help/develop/how-to/packages/python-functions).

Custom model functionality is implemented as a two-step process using train and apply functions. Train function is used
to build/train a model based on provided features and a target variable. The trained model is then stored as an object
inside a given directory. Apply function receives previously saved model and applies it to provided feature columns. The
resulting predictions are returned in a form of a single column in a dataframe.

Custom models utilize [predictive model](learn.md) interface alongside other out-of-the-box solutions.
Once implemented, custom models can be chosen from the list of model engines.

## Train

### Train header parameters

* `#meta.mlname: CustomName` – Name of the custom train function
* `#meta.mlrole: train` – Action (train or apply)
* `#description: Custom ML train function for KNN algorithm` – Description of the function
* `#language: python` – Language (python, r, julia)

### Train input data parameters

* `#input: dataframe df` – Dataframe for training
* `#input: string predictColumn` – List of features/column names separated by a comma

### Train input model parameters

* `#input: int parm1 {category: group1}` - Parameter 1 to build a model. Category is used for input parameters grouping
  within the UI

### Train output model parameters

* `#output: blob model` – Trained model object name. Complete trained model with an absolute path address is stored
  as a blob object to be retrieved and applied by the Apply function

## Apply

### Apply header parameters

* `#meta.mlname: CustomName` – Name of the custom apply function (should match corresponding train function name)
* `#meta.mlrole: apply` – Action (train or apply)
* `#description: Custom ML apply function for KNN algorithm` – Description of the function
* `#language: python` – Language (python, r, julia)

### Apply input model parameters

* `#input: blob model` – Trained model object name. Complete trained model with an absolute path address saved by the
  Train function.

### Apply input data parameters

* `#input: dataframe df` - dataframe for prediction (contains only prediction features)
* `#input: string nameskeys` – optional list of original features/column names separated by comma
* `#input: string namesvalues` – optional list of new features/column names separated by comma. If both `nameskeys`
  and `namesvalues` are supplied, `nameskeys` will replace corresponding `namesvalues` feature names before accessing the
  dataframe

### Apply output parameters

* `#output: dataframe data_out` – single-column dataframe of predicted values

## IsApplicable

`IsApplicable` is a predicate used to determine if a model can be trained to predict the target column using the feature columns. It usually checks the types of input columns as well as their sizes.

Defining `IsApplicable` is required. If it is not defined, the model will not appear in the list of available model engines.

### IsApplicable Header Parameters

* `#meta.mlname: CustomName` – Name of the custom IsApplicable function (should match the corresponding train function name)
* `#meta.mlrole: isApplicable` – Predicate name (isApplicable)
* `#description: Custom ML IsApplicable function for KNN algorithm` – Description of the function
* `#language: python` – Language (python, r, julia)

### IsApplicable Input Data Parameters

* `#input: dataframe df` – Dataframe for training
* `#input: string predict_column` – List of features/column names separated by a comma

### IsApplicable Output Parameters

* `#output: bool result` – Boolean value indicating if the model is applicable to the given data

## IsInteractive

`IsInteractive` is an optional predicate that operates in the same manner as `IsApplicable` but is used to determine if a model can be trained quickly, allowing interactive training.

## Parameter naming convention

**Important:**
All parameter names must **exactly match** those used in the examples.

Do **not** rename, abbreviate, or reformat parameter names — even minor changes can cause errors or unintended behavior.

## Example

### `Train` function

```python
#name: PyKNNTrain
#meta.mlname: PyKNN
#meta.mlrole: train
#description: Custom Python train function for KNN
#language: python
#input: dataframe df
#input: string predictColumn
#input: int n_neighbors {category: FirstParm}
#input: string weights=uniform {category: Parameters; choices: ["uniform", "distance"]}
#input: int leaf_size=30 {category: Parameters}
#input: int p=1 {category: Parameters; range:1-2}
#input: string metric = minkowski {category: Parameters; choices: ["euclidean", "manhattan", "chebyshev", "minkowski"]}
#input: string algorithm = auto {category: Parameters; choices: ["auto","ball_tree", "kd_tree", "brute"]}
#output: blob model

# Import necessary packages
import numpy as np
import pickle
from sklearn.neighbors import KNeighborsClassifier

# Extract/prepare train features and target variable
trainX = df.loc[ :,df.columns != predictColumn]
trainY = np.asarray (df[predictColumn])

# Build and train model
trained_model = KNeighborsClassifier(
    n_neighbors = n_neighbors,
    weights = weights,
    leaf_size= leaf_size,
    p = p,
    metric = metric,
    algorithm = algorithm
)
trained_model.fit(trainX, trainY)

# Save trained model
model = pickle.dumps(trained_model)
```

### `Apply` function

```python
#name: PyKNNApply
#meta.mlname: PyKNN
#meta.mlrole: apply
#description: Custom Python apply function for KNN
#language: python
#input: blob model
#input: dataframe df
#input: string nameskeys [Original features' names]
#input: string namesvalues [New features' names]
#output: dataframe data_out

# Load necessary packages
import numpy as np
import pickle

# If original(nameskeys) and new(namesvalues) passed, map original names to new
nameskeys = nameskeys.split(",")
namesvalues = namesvalues.split(",")
if len(nameskeys) > 0:
    featuresNames = list(df)
    for i in range(len(nameskeys)):
        df = df.rename(columns = {namesvalues[i]: nameskeys[i]})

testX = np.asarray(df)

# Retrieve saved/trained model
trained_model = pickle.loads(model)

# Predict using trained model
predY = trained_model.predict(testX)
data_out = pd.DataFrame({'pred': predY})
```

### `IsApplicable` function


```python

#name: PyKNNIsApplicable
#meta.mlname: PyKNN
#meta.mlrole: isApplicable
#description: Custom Python isApplicable function for KNN
#language: python
#input: dataframe df
#input: string predict_column
#output: bool result


# checks all columns are numerical
import numpy as np
numeric = df.select_dtypes(include=np.number).columns.tolist()
result = len(numeric) == df.shape[1]

```


See also:

* "Custom machine learning models" [video](https://www.youtube.com/watch?v=G66MN30ZPGQ)
