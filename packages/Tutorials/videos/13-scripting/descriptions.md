# Scripting and automation

Playlist: **Scripting and automation**. Feature name of every clip: `How-to · Scripting`. Features in
`public/packages/Tutorials/bdd/features/howto/13-scripting/`.

## Plan

For an analyst who wants to repeat or automate what they do in the UI. JavaScript only: a JavaScript script
runs in the browser; Python, R, Julia and Jupyter run on the server and stay out.

1. The console: open it from the status bar; a dialog's action is logged as a command; the same command typed
   with other values; a command's name alone opens a form for its inputs (demog).
2. Write a JavaScript script: Browse > Platform > Functions > Scripts > New > JavaScript Script; the header
   lines become the run form (inputs with captions and defaults, an output); Run; the result under the editor;
   Save, and the script is found in Scripts.
3. Run a script on your table: a script with a table and two column inputs that adds a column; Run picks the
   open table and the columns; the new column shows in the table (demog).
4. Run a saved script without opening it: its card in Scripts, the context panel's Details, Script (its code)
   and Run panes, the result in the Variables panel.
5. Learn from the sample scripts: Browse > Platform > Functions > JS API Samples, a sample changed and run.
6. Give a script a friendly form: header options (numerical column picker, min/max slider, choices), a script
   that selects the rows far from a column's average (demog).

Not made: "Find and run any function" (Browse > Platform > Functions lists every function, including the
stand's test functions and internal commands, and its search is loose, so a clip would show mostly clutter);
"Get the code for a viewer" (To Script > To JavaScript is a DevTools item whose code needs a table view to
run, which a script run from the editor does not have); "Variables" on its own (the console does not take
column expressions, so there is little to show beyond clip 04's result).

## 01 - The console (0:47)

Title: The Datagrok console: repeat what you did as a command

Description: The console opens from the Console toggle at the bottom right of the status bar. What you do in
the UI shows up there as a command: here Add New Column adds a BMI column and the console logs the
AddNewColumn command with its formula. Typing the same command with another formula and name adds a second
column without the dialog, and typing a command's name alone opens a form for its inputs. Uses the demog demo
dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/01-the-console.feature`

## 02 - Write a JavaScript script (1:10)

Title: Write a JavaScript script in Datagrok

Description: How to write a small script that runs in your browser. Scripts are in Browse under Platform >
Functions > Scripts; New > JavaScript Script opens a template. The comment lines at the top give the script its
name, its inputs (with captions and default values) and its output, and Run turns the inputs into a form. The
example works out a body mass index from height and weight, shows the result under the code, and is saved so
it appears in the Scripts list.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/02-write-a-javascript-script.feature`

## 03 - Run a script on your table (1:04)

Title: Run a JavaScript script on your table in Datagrok

Description: How to run a saved script on a table you have open. The script is found in Browse > Platform >
Functions > Scripts by a search and opened with Edit... from its right-click menu: a dataframe input takes a
table, and a column input takes one of its columns. Run asks for the table and the columns in a form, and the
script adds a BMI column to the demog demo dataset.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/03-run-a-script-on-your-table.feature`

## 04 - Run a saved script without opening it (0:41)

Title: Run a saved script from the context panel in Datagrok

Description: How to read and run a script you saved earlier without opening it in the editor. Find it in
Browse > Platform > Functions > Scripts and click it: the context panel shows what it does, what it takes and
returns, and its code, and its Run pane is the same form you get when you run it from the editor. The result of
the run shows in the Variables panel.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/04-run-a-saved-script.feature`

## 05 - Learn from the sample scripts (0:46)

Title: Learn JavaScript scripting from the Datagrok samples

Description: Datagrok comes with short JavaScript samples for common tasks, in Browse under Platform >
Functions > JS API Samples, grouped by topic. A sample opens as a script you can change and run: here the
color-coding sample, with its height colors changed, opens a demo table with colored height, age, weight and
race columns. Run a script with the play icon above the code or F5.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/05-learn-from-the-samples.feature`

## 06 - Give a script a friendly form (1:13)

Title: Give your Datagrok script a friendly form: sliders, lists, columns

Description: The options in braces after a script's input decide how its run form looks: {type: numerical}
offers only number columns, min and max add a slider, and choices make a drop-down list. The example is a saved
JavaScript script that selects the rows where a value is far from the column's average; run on the demog demo
dataset with WEIGHT and 2.5 standard deviations, it selects 114 rows.

Feature: `public/packages/Tutorials/bdd/features/howto/13-scripting/06-give-a-script-a-friendly-form.feature`
