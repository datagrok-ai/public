# Validation

Validation is performed on study import using [CDISC Rules Engine](https://github.com/cdisc-org/cdisc-rules-engine) against the SENDIG 3.1 standard.

The view consists of two linked tables:
- **Issue Summary** (top): list of violated rules with counts and affected datasets
- **Issue Details** (bottom): detailed violations filtered by the currently selected rule

Selecting a rule in the summary table filters the details table to show only the violations for that rule.

Cells with validation errors are marked with a red indicator. Hover over the mark to see a tooltip with rule ID, message, and affected variables.

For some rules, automatic fixes are available:
1. Select the corresponding rule and click the fix action in the 'action' column
2. Review the proposed fixes in the context panel
3. Click 'Apply fixes' on the ribbon panel to apply the changes

Available fix types include ISO8601 date format corrections and numeric value extraction from string results.

![Validation](validation_view.gif)
