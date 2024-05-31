<!-- TITLE: Tests: Predictive Models Browser -->
<!-- SUBTITLE: -->

# Tests: Predictive Models Browser

Model browser allows you to view previously created [models](learn.md), view their properties and manage
them.

## Testing scenarios

1. Open "Browse Models" from **ML | Models**

* Model browser is open
* You can change the view, use sort and search

1. Use context menu for apply [model](learn.md) from browser

* Models are applied successfully and correctly from context menu

1. Edit [model](learn.md) from its context menu in browser

* Model edit dialog has been opened

1. Share [model](learn.md) from its context menu to another [user](../govern/access-control/users-and-groups#users)
   or user [group](../govern/access-control/users-and-groups#groups)

* Access to the model is enabled for the selected [user](../govern/access-control/users-and-groups#users)
  or [group](../govern/access-control/users-and-groups#groups) of users.
* Appropriate [user](../govern/access-control/users-and-groups#users) received a notification and a letter to the mail

1. Open the "General" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "General" tab is open
* The correct and actually information for all fields is displayed (Created, Created by, Trained on). Structure of the
  model is also shown here

1. Open the "Details" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Details" tab is open
* Display inputs, outputs and using method.
* If a suitable [table](../datagrok/table.md) is opened, the field "Applicable to" is shown

1. Open the "Performance" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Performance" tab is open
* Indicators of the quality of the model shown here (mse, rmse, r2, auc, etc.)
* For OpenCPU models, quality graphs are displayed

1. Open the "Activity" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Activity" tab is open
* Display information about actual actions with [model](learn.md)

1. Select two models (using ```CTRL``` key)

* Selected models appeared on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

1. Click on *"Compare"* action from "Commands" tab on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* New table view "Compare models" added
* Table shows models in rows, and their quality indicators in columns

1. Repeat 9-10 steps with different number of selected models and trained on different servers (
   OpenCPU and H2O)

See also:

* [Predictive modeling](learn.md)
* [Predictive modeling test](predictive-modeling-test.md)
* [Predictive modeling tutorial](../_internal/tutorials/learn.md)
