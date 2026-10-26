# PmxStudio changelog

## 0.0.1 (2026-10-26)

PMX Studio Preclinical, Plan and Analyze, on Crux (crate tgi) and the atlas of historical controls:

* Plan (`pmxPlanStudy`, app route `/plan`): de novo design of a xenograft or syngeneic efficacy study: end of study,
  kill rate per arm, analyzable animals per arm for the target power (assurance over the atlas's labs) with the
  multiplicity correction, animals to randomize and implant, power by group size, end of study and measurement
  schedule, the study calendar and the protocol. Plans saved as `.pmxplan` files (Save plan, Download; opened in Plan
  from Files).
* Analyze (`pmxAnalyzeStudy`, app route `/analyze`): the study analyzed as Plan pre-specifies it (per-animal slope of
  log volume, each treated group against the control), TGI and T/C, each treated animal against its untreated twin,
  the study against its plan, and tolerability from body weight (weight-loss grades, removals, NCI limits).
* The predicted and observed growth and the body weight are drawn by the package's own band chart (`band-chart.ts`,
  echarts): the area between a lower and an upper column, a line through it, one band per group.
* `pmxPriors`, `pmxMapPrior`: the atlas's control-growth priors and the MAP prior through Crux.
