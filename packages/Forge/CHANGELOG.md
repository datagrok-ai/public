# Forge changelog

## v.next

* GROK-20932: Added the EMS schema (models, training runs, applications), method discovery and the Forge app listing methods and models
* GROK-20932: Added training with XGBoost, cross-validation metrics, the Predictive model view, saving models without their training data, and catalog refresh and delete
* GROK-20932: Added applying saved models: Apply dialog, Forge:applyModel, missing-value handling (skip or impute) at training and application, prediction column tag and application records
* GROK-20932: Switched models and training runs to eager promotion (schema 0.1.9): a saved model is a platform entity its author owns
* GROK-20932: Added the model handler with Details, Performance, Activity, Sharing and History panels (the title with the star and the commands; the metrics, activity and methods as grids with tooltips; a copyable seed), model card, tags, the Apply... and Download commands, the Predicted by column panel, Applicable to as a table input with a preset table for Apply, and Compare as forms in the context panel and the Compare view
* GROK-20932: Added every EDA method to the Train view with the suggested default, live retraining (a right-aligned Train button with a tooltip for a method too slow to retrain on every change), data storage modes (reference, copy), a resizable splitter between the inputs and Results, and the `Fix the settings.` hint for an invalid setting
* GROK-20932: Added the Preparation section: one-hot encoding with recorded categories (on by default for up to 20 categories), skipping unique categories (recorded columns), Predict probability with a validated cutoff and AUC-ROC

## 0.0.1 (2026-09-30)
