<!-- TITLE: Console -->
<!-- SUBTITLE: -->

# Console

Use console to call [functions](functions/function.md) 
and record [macros](functions/function.md#macros).

To open: **View | Console**, or press "~" (tilde)

## Controls

|             |              |
|-------------|--------------|
| "~" (tilde) | Open console |
| Tab         | Complete command |
| Up/Down     | Previous/next command |

Two icons on top let you clear the console, or open [variables view](variables-view.md).
Clicking on the function name will bring up its details in the [property panel](property-panel.md).

![Autocomplete](../uploads/gifs/console-autocomplete.gif "Console autocomplete")

## Command examples

```
Mul(2,3)
```
Run `Mul command (multiply two numbers) with the specified parameters.

```
Mul
```
Edit parameters of the `Mul` command and evaluate it in a dialog window.

```
Mul?
```
Get help for the KNN command

```
SelectRows("demog", IsNull("HEIGHT"))
```
Select rows with empty values in the "HEIGHT" column

```
ExtractRows("demog", IsNull("HEIGHT"))
```
Extract rows with empty values in the "HEIGHT" column into a new data frame

## Macros

Whenever a function is executed, it gets logged in the console. Simply copy-and-paste it to execute it again. This can also be used in data transformations and data pipelines.

![](console-macros.gif "Console macros")
 

See also:

  * [Scripting](../develop/scripting.md)
  * [Variables view](variables-view.md)
  * [Grok scripting](grok-script.md)
  * [Scripting plugin](../develop/scripting.md)
