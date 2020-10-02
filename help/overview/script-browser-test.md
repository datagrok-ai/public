<!-- TITLE: Tests: Script Browser -->
<!-- SUBTITLE: -->

# Tests: Script Browser

The script browser allows you to view previously created [scripts](../develop/scripting.md), 
view their properties and manage them.

## Testing scenarios

1. Open "Browse scripts" from (**Tools | Scripting**)
   * Script browser is open
   * You can change the view, use sort and search

1. Run "ACF" [R Script](../develop/scripting.md) from its context menu (input table "TSLA.csv" must be opened)
   * Script runs successfully after selecting the required input parameters

1. Choose "Edit" item from context menu for "ACF" [R Script](../develop/scripting.md)  
   * Script edit tab has been opened

1. Share "ACF" [R Script](../develop/scripting.md) from its context menu to another user or user group
   * Access to the script is enabled for the selected user or group of users
   * Appropriate users received a notification and a letter to the mail

1. Open the "General" tab in [Property Panel](../overview/property-panel.md)
   * "General" tab is open 
   * The correct and actually information for all fields is displayed (Created, Created by, Language)

1. Open the "Details" tab in [Property Panel](../overview/property-panel.md)
   * "Details" tab is open
   * Display inputs and outputs values

1.Open the "Script" tab in [Property Panel](../overview/property-panel.md)
   * "Script" tab is open. The text of the script is shown here. 

1. Open the "History" tab in [Property Panel](../overview/property-panel.md)
   * "History" tab is open
   * Display actually information about running of script
   * Here you can see the status and start time. You can open the [Property Panel](../overview/property-panel.md) for a separate script run here

1. Open the "Statistics" tab in [Property Panel](../overview/property-panel.md)
   * "Statistics" tab is open
   * Display information about runs count, first and last runs

1. Open the "Activity" tab in [Property Panel](../overview/property-panel.md)
   * "Activity" tab is open
   * Display information about actual actions with script

1. Open last run of "ACF" R Script from "History" tab
   * The properties panel for the last run of the script was opened.
   * Here information about the user who made the run, start and end time, duration, input and output values

See also:
  * [Scripting](../develop/scripting.md)
  * [Scripting tutorial](../_internal/tutorials/scripting.md)
  * [Scripting test](../develop/scripting-test.md)
