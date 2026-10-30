# EDA data features — what is still missing

One section per group of a round: the scenarios that need a step the library lacks, the questions for the
help or the product, and the findings.

## EDA dialogs

Feature: `impute-missing-values.feature` (30 Oct 2026, local stand).

### Questions (EDA dialogs)

- **The name of the filled copy.** With In-place off the result is a new column named "b(copy)" (then
  "b(copy 1)"…, `knn-imputer.ts`); the help (`explore/missing-values-imputation.md`) says only "create a new
  column with imputed results". The feature claims the measured name. Should the help give it?
- **Keep empty** shows only when some missing value cannot be imputed (no neighbour has the features), and
  the per-feature weights hide behind the gear next to Distance. The help mentions neither. Not claimed; the
  imputation matrix belongs to EDA `src/tests/mis-vals-imputation-tests.ts`.
- **Nothing to impute**: with every Impute column unticked RUN is disabled and the other inputs disappear
  (`hideWidgets`), with no message. Should the dialog say why? The feature claims RUN disabled only.
