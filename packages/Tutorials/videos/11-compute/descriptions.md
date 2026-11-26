# Scientific computing

Playlist: **Scientific computing**

Diff Studio, Datagrok's app for models made of ordinary differential equations: everything in these clips is
solved, fitted and analyzed in the browser, with no scripts on the server. The audience is a scientist or engineer
who models processes.

## Plan

1. Open and run a model in Diff Studio: Browse > Apps, the hub, a library model, its inputs and the solution
   following them. Done (01).
2. Write your own model: the equations editor, a small model of your own, run it. Done (02).
3. Fit a model to data (parameter optimization): the target data, the parameters to fit, the result. Done (03).
4. Sensitivity analysis: which inputs matter, the plots it gives. Done (04).
5. Find and run models in Model Hub. Done (05).
6. Save a model to your files and open it again. Done (06).
7. Extra: export a model's equations as LaTeX for a paper or a report. Done (07).

Left out: saving to the shared Library (writes into System:AppData/DiffStudio for everyone), sharing a model file
(a grant on the server), Diff Studio's `</>` script export (a script view is its own topic).

## 01 - Run a model in Diff Studio (0:45)

Title: Run a differential equation model in Datagrok Diff Studio

Description: Diff Studio solves ordinary differential equations in the browser. This video opens it from
Browse > Apps, picks the Lotka-Volterra predator-prey model from its library, and changes the starting number of
prey, the prey growth rate (with its slider) and the end time; the curves and the table of values are recomputed
each time. Pointing at an input shows what it stands for.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/01-run-a-model-in-diff-studio.feature

## 02 - Write your own model (1:05)

Title: Write your own ODE model in Datagrok Diff Studio

Description: A Diff Studio model is plain text: a name, the equations, the time span, the starting values and the
parameters. This video starts from the Create button (Browse > Apps > Diff Studio), turns the template into a
logistic population growth model, dN/dt = r * N * (1 - N / K), and solves it with Refresh. With Edit turned off,
every starting value and parameter has become an input with a slider, and dragging the growth rate changes the
curve at once.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/02-write-your-own-model.feature

## 03 - Fit a model to data (1:22)

Title: Fit a model to data in Datagrok Diff Studio

Description: Fitting finds the parameter values that make a model match your data. This video opens the
two-compartment PK model from Diff Studio's library and clicks Fit. The target is a table of drug concentrations
over three doses, picked from Files > App Data > DiffStudio > library (pk.csv, a demo file); the model is set to
three doses as well, and the five PK parameters are switched on to be fitted within their ranges. Run fits in the
browser and lists the fits best first, each with its parameter values and the model drawn over the data.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/03-fit-a-model-to-data.feature

## 04 - Find which inputs matter with sensitivity analysis (1:12)

Title: Sensitivity analysis of a model in Datagrok Diff Studio

Description: Sensitivity analysis varies a model's inputs together and shows which of them change the result the
most. This video opens the two-compartment PK model from Diff Studio's library, clicks Sensitivity, picks the Sobol
method with 1000 samples and switches on the five PK parameters; the output is the drug level 12 hours after the
dose. After Run, the Sobol bar charts show that clearance matters most for the central concentration, and the
correlation plot shows the direction: the higher the clearance, the lower the level.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/04-sensitivity-analysis.feature

## 05 - Run a model from Model Hub (0:40)

Title: Find and run a model in Datagrok Model Hub

Description: Model Hub lists the models you can run on your Datagrok server. This video opens it from
Browse > Apps > Model Hub, points at a card to read what the model does, and runs Ball Flight Simulation with a
double-click. Changing the throw's angle and velocity runs the model again: the trajectory, the maximum distance and
the maximum height update each time. The model's ribbon also has fitting and sensitivity analysis.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/05-run-a-model-from-model-hub.feature

## 06 - Save a model to your files (0:40)

Title: Save a Diff Studio model to your files in Datagrok

Description: This video opens the Lotka-Volterra model from Diff Studio's library (Browse > Apps > Diff Studio)
and saves a copy with the SAVE button on the ribbon, under a name of your choice. The model becomes an .ivp file in
Browse > Files > My files, and clicking the file opens the model in Diff Studio again, with its inputs and charts.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/06-save-a-model-to-your-files.feature

## 07 - Export a model's equations to LaTeX (0:46)

Title: Export model equations to LaTeX in Datagrok Diff Studio

Description: Diff Studio writes a model's equations as LaTeX or Markdown. This video opens the Lotka-Volterra model
(Browse > Apps > Diff Studio), clicks the download icon on the ribbon and picks As LaTeX. The dialog shows the LaTeX
source as you change the options; with Standalone off you get just the part to paste into your own document, which
the copy icon puts on the clipboard. Download saves it as a .tex file instead.

Feature: public/packages/Tutorials/bdd/features/howto/11-compute/07-export-equations-to-latex.feature
