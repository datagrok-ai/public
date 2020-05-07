<!-- TITLE: Tests: Predictive Models Browser -->
<!-- SUBTITLE: -->

# Tests: Predictive Models Browser

Model browser allows you to view previously created [models](predictive-modeling.md), view their properties and manage them.

## Testing scenarios

1. Open "Browse Models" from **Tools | Predictive Modeling**
   * Model browser is open
   * You can change the view, use sort and search

1. Use context menu for apply [model](predictive-modeling.md) from browser
   * Models are applied successfully and correctly from context menu

1. Edit [model](predictive-modeling.md) from its context menu in broser
   * Model edit dialog has been opened
   
1. Share [model](predictive-modeling.md) from its context menu to another [user](../govern/user.md) or user [group](../govern/group.md)
   * Access to the model is enabled for the selected [user](../govern/user.md) or [group](../govern/group.md) of users. 
   * Appropriate [user](../govern/user.md) received a notification and a letter to the mail

1. Open the "General" tab in [Property Panel](../overview/property-panel.md)
   * "General" tab is open
   * The correct and actually information for all fields is displayed (Created, Created by, Trained on). Structure of the model is also shown here

1. Open the "Details" tab in [Property Panel](../overview/property-panel.md)
   * "Details" tab is open
   * Display inputs, outputs and using method. 
   * If a suitable [table](../overview/table.md) is opened, the field "Applicable to" is shown

1. Open the "Performance" tab in [Property Panel](../overview/property-panel.md)
   * "Performance" tab is open
   * Indicators of the quality of the model shown here (mse, rmse, r2, auc, etc.)
   * For OpenCPU models, quality graphs are displayed

1. Open the "Activity" tab in [Property Panel](../overview/property-panel.md)
   * "Activity" tab is open
   * Display information about actual actions with [model](predictive-modeling.md)

1. Select two models (using ```CTRL``` key)
   * Selected models appeared on [Property Panel](../overview/property-panel.md)
   
1. Click on *"Compare"* action from "Commands" tab on [Property Panel](../overview/property-panel.md)
   * New table view "Compare models" added 
   * Table shows models in rows, and their quality indicators in columns

1. Repeat 9-10 steps with different number of selected models and trained on different servers (OpenCPU and H2O)  


See also: 
  * [Predictive modeling](predictive-modeling.md)
  * [Predictive modeling test](../tests/predictive-models-test.md)
  * [Predictive modeling tutorial](../_internal/tutorials/predictive-modeling.md)
  