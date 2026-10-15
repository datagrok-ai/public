# Forge changelog

## v.next

* GROK-20932: Added the EMS schema (models, training runs, applications), method discovery and the Forge app listing methods and models
* GROK-20932: Added training with XGBoost, cross-validation metrics, the Predictive model view, saving models without their training data, and catalog refresh and delete
* GROK-20932: Added applying saved models: Apply dialog, Forge:applyModel, missing-value handling (skip or impute) at training and application, prediction column tag and application records
* GROK-20932: Switched models and training runs to eager promotion (schema 0.1.9): a saved model is a platform entity its author owns
* GROK-20932: Added the model handler with Details, Performance, Activity, Sharing and History panels (the title with the star and the commands; the metrics, activity and methods as grids with tooltips; a copyable seed), model card, tags, the Apply... and Download commands, the Predicted by column panel, Applicable to as a table input with a preset table for Apply, and Compare as forms in the context panel and the Compare view
* GROK-20932: Added every EDA method to the Train view with the suggested default, live retraining (a right-aligned Train button with a tooltip for a method too slow to retrain on every change), data storage modes (reference, copy), a resizable splitter between the inputs and the results, and the `Fix the settings.` hint for an invalid setting
* GROK-20932: Added the Preparation section: one-hot encoding with recorded categories (on by default for up to 20 categories), skipping unique categories (recorded columns), Predict probability with a validated cutoff and AUC-ROC
* GROK-20932: Added Insights & Tips, the results in resizable rows of platform viewers (scatter plots, residuals, distributions, statistics, a ROC curve with the cutoff marked, confusion matrix, wrong predictions, class balance, PC and correlation plots, the method's own viewers; kept while a retrain keeps the data; a tooltip on every panel title; short insights with the numbers and a pointer to their panel in the tooltip, shown as that panel's viewer description) with the metrics in the Performance panel, and Run Evaluation with an Evaluation column, metric values below a threshold or overfitting marked in the metrics grid, a model card of one height and a rich model tooltip
* GROK-20932: Added the table pane Forge models in the built-in Models pane's look (the fitting models with the run icon, the model's tooltip and menu, and Train)
* GROK-20932: Added the catalog's description for the AI assistant
* GROK-20932: Added the model's tooltip to the catalog grid's Name cells
* GROK-20932: Added the Runs grid of the Train view: every training of the view that a change did not stop, newest first; click a run to set its inputs again when they differ from the form's
* GROK-20932: Made model elements draggable onto a Space (the panel title's name, the table pane)
* GROK-20932: Added the Tags filter to the Domains view of models (schema 0.1.13)
* GROK-20932: Updated the schema descriptions: training run columns match the model's, and Data source mentions the OpenFile reference
* GROK-20932: Fixed the Train view keeping a column removed from the table in Target or Features
* GROK-20932: Added the API functions Forge:trainModel (the Train view's checks and defaults, every training recorded, optional save with a storage mode) and Forge:deleteModel, and the training run's source (schema 0.1.14)
* GROK-20932: Added readiness attributes for automation: aria-busy on the Train and Evaluation results; aria-invalid and data-forge-run-id (the training run shown) on the Train results; aria-expanded and aria-invalid on the collapsible groups
* GROK-20932: Added the BDD features of the Train view (bdd/): the iris and winequality journeys, each training claimed by the run it shows and that run on the server

## 0.0.1 (2026-09-30)
