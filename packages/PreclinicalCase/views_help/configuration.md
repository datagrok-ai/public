# Configuration

Treatment and control arm mapping for comparative analysis.

This view allows you to define treatment-to-control arm pairings. Arm codes and names are extracted from the TA (Trial Arms) domain.

Once treatment/control pairs are configured, additional computed columns are added to the LB domain:
- **CONTROL_MEAN**: mean of control group values per test and visit day
- **DELTA_VS_CONTROL**: difference between the current value and control mean
- **PCT_VS_CONTROL**: percent difference from control mean

These columns enable direct comparison of treatment effects relative to the control group.

![Configuration](configuration_view.png)
