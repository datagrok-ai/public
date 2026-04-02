# Measurements

Comprehensive analysis of all measurement findings across SEND domains.

This view combines data from all findings domains (LB, BW, CL, EG, VS, PC, etc.) into a unified table with standardized columns:
- **test**: test name (with category/specimen when available)
- **result**: numeric result value
- **string_result**: character result
- **units**: result units
- **visit_day**: study day
- **category**: test category
- **specimen**: specimen type
- **blfl**: baseline flag

Computed columns are added automatically:
- **baseline**: baseline value per subject/test (uses BLFL="Y" flag per SENDIG v3.1.1, falls back to earliest study day)
- **change**: change from baseline (current - baseline)
- **pct_change**: percent change from baseline
- **max_post_value**: maximum value across all post-baseline visits
- **min_pct_change** / **max_pct_change**: extreme percent changes post-baseline

![Measurements](measurements_view.gif)
