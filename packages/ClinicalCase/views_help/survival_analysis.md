# Survival analysis

This view is designed to create Kaplan-Meier curves and perform covariates analysis.

Steps to perform survival analysis:

- create dataset (`Dataset` tab)

Choose endpoint for which you want to create Kaplan-Meier curve. Optionally you can include basic covariates(age, sex, race, treatment arm) into dataset. Click `Create dataset` button.

Dataset can be further filtered.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/survival_dataset.gif" height="500" width='800'/>

- go to `Survival Chart` tab. You will see Kaplan-Meier curve. You can modify confidence interval or choose strata(stratas will be available in dropdown list in case covariates for dataset are selected). After modifying parameters curve will be updated.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/survival_kaplan_meier.gif" height="500" width='800'/>

In case you want to perform covariates analysis you should go to `Covariates` tab and check one or several covariates.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/survival_covariates.PNG" height="500" width='800'/>
