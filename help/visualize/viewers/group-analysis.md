---
title: "Group analysis viewer"
---

Group analysis viewer groups data by one or more columns and analyzes them using aggregations, charts, and statistical tests. It helps compare groups and understand their characteristics.

>Note: To use a group analysis viewer, install the package
[Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts).

## Add a group analysis viewer

1. On the menu ribbon, click the **Add viewer** icon. A dialog opens.
1. In the dialog, select **Group Analysis**.

> Developers: To add the viewer from the console, use:
 `grok.shell.tv.addViewer('Group Analysis');` 

## Configuring a group analysis viewer

* **Select grouping columns**: choose one or more columns under **Group by**. A row is created for each unique combination of values.
* **Add analysis columns**: 
  1. Click the **+** icon to add columns for analysis. A dialog opens.
  1. In the dialog, select:
      * **Column**: The column to analyze
      * **Column type**: `Aggregate`, `Chart`, or `Statistics`
      * **Function/Chart/Statistic**: Specific analysis method
  1. Click OK
  
  The grid updates automatically to show analysis results. Charts appear as interactive visualizations within cells.

   ![Group analysis viewer](img/group-analysis.gif "Group analysis viewer")

## See also
  
- [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts)  
- [Viewers](viewers.md)