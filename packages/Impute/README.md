# Impute
Impute is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platfrom.
It provides the user with a selection of tools designed to impute an incomplete dataset in an orderly predictive fashion, 
so as to retain the rest of the feature values associated with the missing instances. 

# Usage guide
This package is structured as a two-step wizard. Firstly, the user selects which columns they wish to impute and which 
columns they wish to be the source of inference. Upon completing the preliminary visual analysis (assisted by matrix,
cluster and correlation plots) and discarding the undesired columns the user is presented with a choice of algorithms.
These are split into three groups as follows:

* **For numeric data**: applicable to both discrete and continuous variables (but not factors);
* **For categorical data**: applicable to all categorical variables, both string and numeric;
* **For mixed data**: applicable to any combination of variable data types; 

There is no need to preprocess or factorize the variables in advance. All the transformations (along with normalization
and centering, where applicable) are performed internally. If there is a requiremenet for user input it will be explicitly
stated so under the '*algorithm parameters*' tab. Some algorithms might still perform their function if the wrong data type
is selected, however the quality of the resulting inference will be sub-standard. 

# Methods summary
Some of these methods work better than others, depending on the missingness pattern, feature data types and the hidden 
deletion mechanisms, etc. A clear understanding of the dataset and all of its aforementioned properties is required 
to achieve reliable results. Please carefully review the presented selection of models and pay close attention 
to the associated hyper-parameters and their effects on data (note: some adjustable elements are set to sensible 
defaults and hidden from the user).

| Algorithm | Numeric data | Categorical data | Mixed data |
| ---- | ---- | ---- | ---- |
| [Hmisc::aregImpute:pmm](https://www.rdocumentation.org/packages/Hmisc/versions/4.4-1/topics/aregImpute) | &check; | &check; | &check; |
| [Hmisc::aregImpute:regression](https://www.rdocumentation.org/packages/Hmisc/versions/4.4-1/topics/aregImpute) | &check; |  |  |
| [Hmisc::aregImpute:pmm](https://www.rdocumentation.org/packages/Hmisc/versions/4.4-1/topics/aregImpute) | &check; | |  |
| [VIM::kNN](https://www.rdocumentation.org/packages/VIM/versions/6.0.0/topics/kNN) | &check; | &check; | &check; |
| [mice::mice](https://www.rdocumentation.org/packages/mice/versions/3.9.0/topics/mice) | &check; | &check; | &check; |
| [missMDA::PCA](https://www.rdocumentation.org/packages/missMDA/versions/1.17/topics/imputePCA) | &check; | |  |
| [missMDA::FAMD](https://www.rdocumentation.org/packages/missMDA/versions/1.17/topics/imputeFAMD) | &check; | &check; | &check; |
| [missMDA::MCA](https://www.rdocumentation.org/packages/missMDA/versions/1.17/topics/imputeMCA) | | &check; |  |
| [pcaMethods::nipals](https://www.rdocumentation.org/packages/pcaMethods/versions/1.64.0) | &check; | |  |
| [pcaMethods::ppca](https://www.rdocumentation.org/packages/pcaMethods/versions/1.64.0) | &check; | |  |
| [pcaMethods::bpca](https://www.rdocumentation.org/packages/pcaMethods/versions/1.64.0) | &check; | |  |
| [pcaMethods::nlpca](https://www.rdocumentation.org/packages/pcaMethods/versions/1.64.0) | &check; | |  |
| [missForest::missForest](https://www.rdocumentation.org/packages/missForest/versions/1.4) | &check;| &check; | &check; |
