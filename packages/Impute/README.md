# Impute
Impute is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platfrom.
It provides the user with a selection of tools designed to impute an incomplete dataset in an orderly predictive fashion, 
so as to retain the rest of the feature values associated with the missing instances. 

# Usage guide
This package is structured as a two-step wizard. Firstly, the user selects which columns he wishes to impute and which 
columns he wishes to be the source of inference. Upon completing the preliminary visual analysis (assisted by matrix,
cluster and correlation plots) and discarding the undesired columns the user is presented with a choice of algorithms.
These are split into three groups as follows:

* **For continuous data**
* **For categorical data**
* **For mixed data**  

There is no need to preprocess or factorize the variables in advance. All the transformations (along with normalization
and centering, where applicable) are performed internally. If there is a requiremenet for user input it will be explicitly
stated so under the 'algorithm parameters' tab.  

# Methods summary
Some of these methods work better than others, depending on missingness patterns, feature data types and the hidden 
deletion mechanisms, etc. A clear understanding of the dataset and all of its aforementioned properties is required 
to achieve dependable results. Please carefully review the presented selection of models and pay close attention 
to the associated hyper-parameters and their effects on data.

| Algorithm | Continuous data | Categorical data | Mixed data |
| ---- | ---- | ---- | ---- |
| [Hmisc::aregImpute:pmm]() | &check; | &check; | &check; |
| [Hmisc::aregImpute:regression]() | &check; |  |  |
| [Hmisc::aregImpute:pmm]() | &check; | |  |
| [VIM::kNN]() | &check; | &check; | &check; |
| [mice::mice]() | &check; | &check; | &check; |
| [missMDA::PCA]() | &check; | |  |
| [missMDA::FAMD]() | &check; | &check; | &check; |
| [missMDA::MCA]() | | &check; |  |
| [pcaMethods::nipals]() | &check; | |  |
| [pcaMethods::ppca]() | &check; | |  |
| [pcaMethods::bpca]() | &check; | |  |
| [pcaMethods::nlpca]() | &check; | |  |
| [missForest::missForest]() | &check;| &check; | &check; |
