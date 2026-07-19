# Plates changelog

## v.next

* Analyses: break circular import between `base-analysis.ts` and `package.ts` to fix TDZ error `Cannot access 'AnalysisBase' before initialization` when the test bundle loads.
