<!-- TITLE: Tests: Scripting -->
<!-- SUBTITLE: -->

# Tests: Scripting

[Scripting](scripting.md) is an integration mechanism with R, a language for statistical computing.

R Scripting combines fast interactive visualizations and other features of the system with thousands of statistical
packages and visualizations available in R.

The system allows to create your own scripts, open scripts saved locally, and also use ready-made samples

## Testing scenarios

1. Open "New Script..." from **Tools | Scripting**

* The script creation tab is open

1. Click on the "Open template script" icon

* Script creation template is open
* The template has fields "#name", "#description", "#help", "#language", "#tags", "#input", "#output"

1. Click on the "Clear script" icon

* Script input area has become empty

1. Select sample "ACF" from R scripts samples

* Sample "ACF" was opened

1. Click on the "Load sample table" icon

* Table with the data for the sample was opened

1. Click on the "Run script" icon

* Dialog for selecting input script parameters was opened
* Current history of input parameters is present here

1. Run script with incorrect input data types

* Notification of incorrect input data type

1. Run script with correct input data types

* Script runs successfully. The output is correct

1. Remove optional fields from the script ("#description", "#help", #tags") and run it

* Script runs successfully without optional fields

1. Close all open tables and run the script

* Script does not run. Notification about necessary table

1. Run JS script "PCA"

1. Run Python script "Butina Molecules Clustering" on Jupyter server

1. Run Julia script "PCA" on Jupyter server

1. Run Grok script "Demo Scripts Run Test"

Use this test scenario to run scripts on OpenCPU and Jupyter servers.

See also:

* [Scripting](scripting.md)
* [Scripting tutorial](../_internal/tutorials/scripting.md)
* [Script browser test](../overview/script-browser-test.md)
* [Scripting auto test](scripting-test.side)
