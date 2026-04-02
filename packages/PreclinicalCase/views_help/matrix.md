# Matrix Table

A comprehensive matrix view that pivots findings data by animal and study day.

Test results from all findings domains are combined and pivoted so that each column represents a specific test, and each row represents a subject at a particular visit day. Values show the average numeric result for each subject/test/day combination.

The view includes:
- **Matrix Plot** viewer showing pairwise relationships between tests
- Color coding by treatment arm (ARMCD) when available (can be changed via viewer settings)
- Demographic columns joined from DM domain for filtering

This view is particularly useful for identifying correlations across multiple tests and detecting patterns across animals simultaneously.

![Matrix view](matrix_view.gif)
