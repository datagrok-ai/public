---
created: 2026-03-10T14:16:05.170Z
title: Challenge imputation to leverage preexisting DG methods and UI
area: analysis
files:
  - packages/Proteomics/src/analysis/imputation.ts
---

## Problem

The Proteomics package implements imputation methods (MinProb, kNN, zero, mean, median) from scratch in TypeScript. Datagrok may already provide built-in imputation methods, missing value handling utilities, or UI patterns (e.g., MissingValuesImputer, column statistics, or library functions in @datagrok-libraries/ml or @datagrok-libraries/statistics) that could be reused. Custom implementations risk diverging from platform conventions and missing optimizations the platform already provides.

## Solution

Audit Datagrok's existing capabilities:
1. Check `@datagrok-libraries/ml` and `@datagrok-libraries/statistics` for imputation functions
2. Check `DG.DataFrame` and `DG.Column` APIs for built-in missing value handling
3. Check if Datagrok has a standard imputation UI pattern or dialog that could be extended
4. For methods where DG equivalents exist, refactor to wrap/call DG functions instead of custom implementations
5. For proteomics-specific methods (MinProb downshift, valid-values filtering), keep custom but ensure they compose well with DG primitives
