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
* GROK-20932: Added catalog readiness attributes for automation: aria-busy (from a change of the models or a Refresh click, before the reload) and data-forge-catalog-load (the number of the load shown) on the models pane, and the forge-app class on the Forge view
* GROK-20932: Added the BDD features of the Apply dialog, the catalog, the app and the menus (bdd/)
* GROK-20932: Fixed the shown model's Details (Last run, Applications) and Activity staying stale after an application until another row was clicked
* GROK-20932: Fixed the Activity and Performance grids of the model's context panel being narrower than the panel
* GROK-20932: Added the BDD features of preparation and training with gaps, the Train view on demog and the table pane Forge models (bdd/)
* GROK-20932: Fixed the Train view's results keeping aria-invalid after a cancelled training was followed by a completed one
* GROK-20932: Fixed a click on a model's name in the table pane Forge models doing nothing right after the pane was expanded
* GROK-20932: Added the API functions Forge:listMethods, listModels, getModel, updateModel, evaluateModel and shareModel
* GROK-20932: Fixed the model's context panel not following an edit or delete made outside the catalog
* GROK-20932: Added the package knowledge for Datagrok AI (agents/)
* GROK-20932: Added the missing-values mode keep for methods that handle missing values themselves (meta.mlnativenans): XGBoost, Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost and Decision Tree from EDA always train with keep unless skip or impute is chosen explicitly (API), and Missing values is hidden for them in the Train view; the Performance line "Missing values: kept (the method handles them)" and the Runs cell "keep missing values" show only when the training features had gaps
* GROK-20932: Made Apply, Forge:applyModel, Run Evaluation and Forge:evaluateModel handle missing values in the mode the model was trained with (skip, impute with its settings, or keep); a model without a recorded mode keeps the gaps when its method handles missing values and skips those rows otherwise; the Apply dialog presets Missing values to that mode and hides it for keep
* GROK-20932: Added AUC-ROC for every classifier whose method returns class probabilities (all EDA classifiers but SVM; proba_<class> columns of apply; the macro one-vs-rest average for more than two classes) and the ROC Curve for two-class classifiers without Predict probability, in training, Run Evaluation and Forge:evaluateModel; Predict probability stays for the cutoff (schema 0.1.15)
* GROK-20932: Added Balanced accuracy and MCC to the metrics of every classifier and Average precision from Crux (EDA's classificationReport), per-class values saved with the model (schema 0.1.16) and shown in the Per-class panel for three or more classes, the ROC Curve drawn from Crux's exact curve and a PR Curve panel for two classes, in training, Run Evaluation and the API
* GROK-20932: Fixed Predict probability switching a two-class training to Linear Regression (the suggestion read the prepared 0 / 1 target) and the empty class the cut labels showed in Confusions
* GROK-20932: Improved the results: the Per-class panel is a parallel coordinates plot in the target's class colors, without in-chart filters or a color selector, its legend on top when short and hidden otherwise; ROC Curve and PR Curve are thin lines with the cutoff circled, axes FPR (1 − Specificity), Recall (Sensitivity) and Precision, the PR Curve's non-informative baseline at the positive class's share; Performance moved to the second row and Wrong predictions to its own last row above Runs; the metric names Recall (Sensitivity) and Negative Predictive Value

## 0.0.1 (2026-09-30)
