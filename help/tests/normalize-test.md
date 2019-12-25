<!-- TITLE: Tests: Normalize -->
<!-- SUBTITLE: -->

# Tests: Normalize

Performs numerical column normalization (min-max, z-scores)

## Testing scenario


1. Open "demog" table 

1. Open "Normalize..." dialog for *"Age"* column (from it's context menu or from "Actions" tab on [Property Panel](../features/property-panel.md))
   * "Normalize" dialog is open

1. Select categorical columns for *Column* field (*"Sex"*, *"Control"*, *"Started"*)
   * Field is highlighted in red
   * "OK" disabled for clicking
   * Reason is shown in tooltip (column is not numerated)

1. Select *"Age"* in *Column* field and leave default method (min-max) and execute dialog
   * *"Age"* column values ​​are normalized by min-max method
   
1. Open "Normalize..." dialog for *"Height"* column (from it's context menu or from "Actions" tab on [Property Panel](../features/property-panel.md))
   * "Normalize" dialog is open   
   
1. Select "z-scores" for *Method* field and execute dialog
   * *"Height"* column values ​​are normalized by z-scores method
   
1. Open **Tools | Console** (or use hotkey "~")
   * "Console"  is open
   * Console shows Normalize() functions performed in previous steps
   
1. Execute command "Normalize("demog", "Weight")" in console
   * *"Weight"* column values ​​are normalized by min-max method
  