<!-- TITLE: Tests: Console -->
<!-- SUBTITLE: -->

# Tests: console

[Console](../navigation/panels/panels.md#console) is used to execute commands and record macros.

## Testing scenarios

1. Open dataset "demog".

1. Open [console](../navigation/panels/panels.md#console) by clicking ```"~"``` (tilde) button on keyboard.

1. Use menu "Select" for select rows ([Random](../../explore/select-random-rows.md)
   , [Duplicates](../../explore/select-duplicates.md),
   [Missing Values](../../explore/missing-values-imputation.md))

* Commands for opening dialogs and selection operations were written to the console.

1. Add viewers.

* Commands for adding all viewers were written to the console.

1. Show different windows from menu "View". ([Tables](../concepts/table.md), [Browse](../navigation/views/browse.md), etc.)

* Commands for opening all windows were written to the console.

1. Use "Data" and "Tools" menus items for operations with data. ([Split column](../../transform/text-to-columns.md),
   [Categorize](../../transform/categorize-data.md), [Anonymize](../../transform/anonymize-data.md)
   , [Cluster](../../explore/cluster-data.md), etc.)

* Commands for opening dialogs and data operations (with parameters) were written to the console.

1. Run [scripts](../../compute/scripting/scripting.mdx) from Script Browser.

* The command of script run was written to the console.
* The command matches the script name and has selected parameters.

1. Start writing the name of the [script](../../compute/scripting/scripting.mdx) "Linear regression" in the input field of the
   command [console](../../datagrok/navigation/panels/panels.md#console). After entering the "Line"
   press the ```"Tab"``` key from the keyboard.

* After pressing the ```"Tab"``` key, the command was written automatically.

1. Enter the parameters ("demog", "HEIGHT", "WEIGHT", false) for the command "Linearregression" and execute it.

* Script "Linear regression" is run with the entered parameters. Execution of the is successful.
* The result is shown in the console.

1. To navigate through the commands executed, use the ```"Up"``` and ```"Down"``` keys on the keyboard.

1. Highlight and drag the command from the console to the [console](../navigation/panels/panels.md#console) input field.

* Selected command is dragged into the console field.

1. Use symbol "?" in front of the command for get help on it. (?Linearregression)

1. Use nested commands to build macros. (SelectRows("demog", IsNull("HEIGHT"))

* Nested commands are working correctly.
* After starting the example, all rows were selected for which the value "HEIGHT" is empty.

1. Click on "Clear" button

* Console cleared

See also:

* [Console](../navigation/panels/panels.md#console)
* [Scripting](../../compute/scripting/scripting.mdx)
