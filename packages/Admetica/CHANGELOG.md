# Admetica changelog

## 1.0.2 (2024-10-18)

### Bug Fixes

* Admetica: Fixed tests
* Admetica: Support CPU usage when no GPU is available

## 1.0.1 (2024-09-30)

### Features

* Admetica: Introduced toxicity models (hERG and LD50 for acute toxicity)

### Bug Fixes

* Admetica: Fixed an issue where the pie chart in the summary pane failed to generate when editing the structure in the sketcher
* Admetica: Added color coding for calculated values in the property panel for improved clarity
* Admetica: Introduced a new input for selecting color coding rules from templates
* Admetica: Added missing color coding rules for Half-Life model
* Admetica: Resolved an issue with classification models (e.g., metabolism models) incorrectly returning absolute zero as a result
* Admetica: Updated the color scheme by switching to a categorical palette for better visual distinction

## 1.0.0 (2024-09-18)

### Features

* Admetica: Ability to calculate ADMET properties for the entire column
* Admetica: Ability to calculate ADMET properties for individual cells
* Admetica: Seamless transition from Chemprop v1 to v2
* Admetica: Update models aligned with the latest Chemprop v2 features


### Bug Fixes

* GROK-16464: Chem | Adme: Form does not work in the demo
* GROK-14645: Chem | Adme: Errors when trying to calculate Adme/Tox on the Context pane
* GROK-13458: Chem | Adme: Incorrect handling of malformed data
* GROK-13456: Chem | Adme: Pgp-substrate wrong color coding
* GROK-13024: Chem | Adme: Multiple errors in the property panel