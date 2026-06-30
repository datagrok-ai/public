---
created: 2026-03-10T14:16:05.170Z
title: Challenge normalization to leverage preexisting DG methods
area: analysis
files:
  - packages/Proteomics/src/analysis/normalization.ts
---

## Problem

The Proteomics package implements normalization methods (median centering, quantile, VSN) from scratch. Datagrok or its libraries (@datagrok-libraries/statistics, @datagrok-libraries/ml) may already provide normalization utilities, column scaling functions, or statistical transforms that could be reused. Custom implementations increase maintenance burden and risk diverging from platform conventions.

## Solution

Audit Datagrok's existing capabilities:
1. Check `@datagrok-libraries/statistics` for normalization or scaling functions
2. Check `@datagrok-libraries/ml` for data preprocessing/normalization utilities
3. Check `DG.Column` API for built-in statistical operations (median, quantile computation)
4. Check if Datagrok has standard normalization patterns used by other packages (e.g., Chem, Bio)
5. For methods where DG equivalents exist, refactor to wrap/call DG functions
6. Keep VSN as R-script (domain-specific, no DG equivalent expected) but verify R script conventions match platform standards
