# Diff Studio

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.09090/status.svg)](https://doi.org/10.21105/joss.09090)
[![Run: Diff Studio](https://img.shields.io/badge/app-%20Diff%20Studio-red.svg)](https://public.datagrok.ai/apps/DiffStudio)
[![Tutorial: Differential equations](https://img.shields.io/badge/tutorial-Differential%20Equations-yellow.svg)](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://datagrok.ai/help/compute/diff-studio)
[![License: MIT](https://img.shields.io/badge/license-MIT-success.svg)](https://github.com/datagrok-ai/diff-grok/blob/main/LICENSE)

Diff Studio is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai) platform. It provides in-browser tools for solving [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem) for systems of ordinary differential equations (ODE) and its interactive exploration.

* Go to **Apps > Compute** and run **Diff Studio**
* Enter formulas or modify template
* Run computations and explore your model

The solver takes a set of the differential equations in a declarative form, and creates a UI that solves the equations, visualizes the results, and lets you change parameters on the fly.

## Quick start

Datagrok provides intuitive tools for the rapid solving ODEs.

* Launch Diff Studio:
  * Go to **Apps > Compute** and run **Diff Studio**
  * The **Hub** opens with ready-to-run models grouped into **Templates**, **Library** and **Recent**
* Run a model:
  * Double-click any card on the Hub to run it
  * Or click **Create** to start from a basic template
  * Or click **Upload** to load an *ivp*-file from your disk
  * Or drag an *ivp*-file right into the browser
  * Modify inputs and explore computations
* Modify equations:
  * Turn on the **Edit** toggle on the top panel
  * Edit formulas or add new ones
  * Click <i class="fas fa-sync"></i> **Refresh** or press **F5** to apply changes
* Save model:
  * Click the **Save** button on the top panel to save model to your platform files (**Browse > Files > My files**)
  * Click <i class="fas fa-arrow-to-bottom"></i> **Download** to save model to a local file. Find the *ivp*-file in Downloads. You can open and edit this file using any text editor
  * Click <i class="fas fa-layer-plus"></i> **Save to Library** to publish the model to your personal **Library** — it appears on the Hub and in the browse tree
* Reopen a model:
  * Pick it from the **Hub** (Templates, Library, Recent) — double-click the card
  * Or use the browse tree: **Browse > Apps > Compute > Diff Studio > Templates | Library | Recent**
  * Or use the <i class="fas fa-folder-open d4-combo-popup" style="min-width: 0px; cursor: default"></i> **Open** menu on the top panel — it also exposes **Import...** and **My Models**
* Analyze model
  * Turn off the **Edit** toggle on the top panel
  * Click the **Fit** icon to [optimize inputs](https://datagrok.ai/help/compute/function-analysis#parameter-optimization)
  * Click the **Sensitivity** icon to run [sensitivity analysis](https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis)

## Hub

The **Hub** is the default landing view of Diff Studio. It shows all available models as cards and exposes the most common actions.

* Action buttons at the top:
  * <i class="fas fa-sync"></i> **Refresh** — reload the Hub contents
  * **Create** — start a new model from the basic template
  * **Upload** — load a local *ivp*-file and run it
* Card sections:
  * **Templates** — starter models (`Basic`, `Advanced`, `Extended`)
  * **Library** — built-in use cases plus models you have saved to your personal Library
  * **Recent** — models you have opened recently, including custom files
* Card interactions:
  * Hover a card for a tooltip with the model description
  * Double-click a card to run the model
  * Right-click a card for a context menu: **Run**, **Copy link**, **Help**. Custom Library cards also expose **Settings...** for editing the help link

The same sections are mirrored in the browse tree under **Browse > Apps > Compute > Diff Studio**.

## Create model from template

Start from one of these templates:

| Template | Features|
|----------|---------|
| `Basic`    | the simplest model|
| `Advanced` | extra math features, including *expressions*, *constants*, *parameters* and *tolerance* specification|
| `Extended` | the *annotating* feature for extended UI generation                 |

## Library

The **Library** section on the Hub contains built-in example models. You can run them directly from Hub cards or from the browse tree, and use them as a starting point — they cover the main features of Diff Studio.

* `Chem reactions`
  * simulates [mass-action kinetics](https://en.wikipedia.org/wiki/Law_of_mass_action)
  * illustrates annotation of inputs
* `Robertson's model`
  * Robertson's chemical reaction model
  * [stiff equations](https://en.wikipedia.org/wiki/Stiff_equation) example
  * shows how Datagrok solves complicated ODEs
* `Fermentation`
  * models the kinetics of the biochemical reactions in [fermentation](https://en.wikipedia.org/wiki/Fermentation)
  * shows the usage of `min` and `max` in inputs annotation
* `PK`
  * [pharmacokinetic](https://en.wikipedia.org/wiki/Pharmacokinetics) (PK) simulation
  * demonstrates the usage of the `meta.solver` feature for numerical solver management
* `PK-PD`
  * simulates pharmacokinetics (PK), pharmacodynamics (PD), and their [relationship](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046)
  * illustrates the usage of the `loop` feature for dosing specification
* `Acid production`
  * models gluconic acid [production](https://oatao.univ-toulouse.fr/9919/1/Elqotbi_9919.pdf) by Aspergillus niger
  * shows the usage of the `update` feature for multi-stage simulation
* `Nimotuzumab`
  * models population pharmacokinetic for [nimotuzumab](https://www.mdpi.com/1999-4923/12/12/1147)
  * demonstrates the `output` feature
* `Bioreactor`
  * simulates the bioreactor [processes](https://doi.org/10.1074/jbc.RA117.000303)
  * shows how to use `meta.inputs` to specify a table with pre-defined model inputs
* `Pollution`
  * describes a chemical reaction part of the air pollution [model](https://archimede.uniba.it/~testset/report/pollu.pdf) consisting of 25 reaction and 20 reacting compounds
  * demonstrates the simulation of processes described by a stiff system of ODEs

### Custom models

You can extend the **Library** with your own models:

* Open or write a model in the editor, then click <i class="fas fa-layer-plus"></i> **Save to Library** on the top panel. The file is stored under `System:AppData/DiffStudio/library` and registered in the Library manifest (`external-models.json`).
* The saved model appears on the Hub in the **Library** section and in **Browse > Apps > Compute > Diff Studio > Library**.
* Right-click a custom model card and choose **Settings...** to set a help link for it.

Datagrok's ODEs suite has tools for solving both [stiff](https://en.wikipedia.org/wiki/Stiff_equation) and non-stiff equations. It provides a [numerical solution](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations). The primary solvers are [CVODE](https://sundials.readthedocs.io/en/latest/cvode/index.html) and [LSODA](https://doi.org/10.1137/0904010) — variable-order methods that automatically detect stiffness and switch between Adams (non-stiff) and BDF (stiff) formulations.

## Model structure

A model defines [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem). It contains *name*, *differential equations*, *initial values* and *argument* specifications:

```python
#name: My ODEs
#equations:
  dx/dt = x + y + exp(t)
  dy/dt = x - y - cos(t)
#argument: t
  initial = 0
  final = 1
  step = 0.01
#inits:
  x = 2
  y = 5
```

Use the following sections to specify various problems:

|Control block|Features|
|-------------|--------|
|```#name```|Defines a name|
|```#equations```|Differential equation specification|
|```#argument```|Independent variable specification|
|```#inits```|Initial values specification|
|```#constants```|Constants specification|
|```#parameters```|Parameters specification|
|```#expressions```|Defines auxiliary compuations|
|```#output```|Defines output columns and their captions|
|```#description```|Defines description of the model|
|```#comment```|Specifies comments block|
|```#meta.solver```|Defines the solver settings|
|```#meta.inputs```|Path to the table with pre-defined model inputs|

## Annotations

The **Run** tab of Diff Studio provides UI for interactive model exploration. Use [annotations](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation#parameter-types-and-options) to improve usability. Append them, when defining *parameters*, *initial values* and *argument*.

Group inputs by specifying their `category`:

```python
#parameters:
  P1 = 1 {category: Parameters}
  P2 = -1 {category: Parameters}
```

Define the desired `caption`:

```python
#argument: t
  start = 0 {caption: Initial; category: Time}
  finish = 2 {caption: Final; category: Time}
  step = 0.01 {caption: Step; category: Time}
```

Specify `min`, `max` and `step` values to get sliders and clickers for the rapid model exploration:

```python
#inits:
  x = 2 {min: 0; max: 5}
  y = 0 {min: -2; max: 2; step: 0.1}
```

Provide hints in brackets `[ ]`:

```python
  P1 = 1 {category: Parameters} [P1 parameter tooltip]
```

## Solver management

You can manage the solver of ODEs. Specify its settings in the `#meta.solver`-line:

* the numerical method (`method`)
* the maximum number of iterations (`maxIterations`)
* the maximum computation time (`maxTimeMs`)

Diff Studio implements the following numerical methods for solving ODEs:

Automatic stiffness-detecting methods:

|Method|Value|
|-------------|--------|
|[CVODE](https://sundials.readthedocs.io/en/latest/cvode/index.html) - variable-order, variable-step BDF solver from [SUNDIALS](https://github.com/LLNL/sundials) ([Hindmarsh et al., 2005](https://doi.org/10.1145/1089014.1089020); [Cohen & Hindmarsh, 1996](https://doi.org/10.1063/1.4822377))|`'cvode'`|
|[LSODA](https://doi.org/10.1137/0904010) - variable-order solver with automatic switching between Adams (non-stiff) and BDF (stiff)|`'lsoda'`|

Implicit methods (for [stiff](https://en.wikipedia.org/wiki/Stiff_equation) ODEs) - [Rosenbrock-Wanner](https://doi.org/10.1016/j.cam.2015.03.010) type:

|Method|Value|
|-------------|--------|
|the modified Rosenbrock triple ([MRT](https://doi.org/10.1016/j.cam.2015.03.010))|`'mrt'`|
|the [ROS3PRw](https://doi.org/10.1016/j.cam.2015.03.010) method|`'ros3prw'`|
|the [ROS34PRw](https://doi.org/10.1016/j.cam.2015.03.010) method|`'ros34prw'`|

Explicit methods (for non-stiff ODEs) - [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) type:

|Method|Value|
|-------------|--------|
|the Bogacki-Shampine 3(2) method ([RK3](https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method))|`'rk3'`|
|the Runge-Kutta-Fehlberg 4(5) method ([RK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method))|`'rk4'`|
|the Dormand-Prince 5(4) method ([RKDP](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method))|`'rkdp'`|

Explicit methods (for non-stiff ODEs) - [Adams-Bashforth](https://en.wikipedia.org/wiki/Linear_multistep_method) type:

|Method|Value|
|-------------|--------|
|the predictor-corrector method of order 4 ([AB4](https://en.wikipedia.org/wiki/Linear_multistep_method))|`'ab4'`|
|the predictor-corrector method of order 5 ([AB5](https://en.wikipedia.org/wiki/Linear_multistep_method))|`'ab5'`|

By default, Diff Studio uses ROS34PRw. You may specify the method as follows:

```python
#meta.solver: {method: 'cvode'}
```

To check correctness of formulas, set the maximum number of iterations:

```python
#meta.solver: {method: 'mrt'; maxIterations: 1}
```

Diff Studio alerts you if computations take too long. The default time limit is 5 seconds. To customize it, set the maximum computation time (in milliseconds):

```python
#meta.solver: {method: 'mrt'; maxTimeMs: 50}
```

## Lookup tables

Lookup tables are pre-defined sets of model input values. They're organized as follows:

||x|y|...|
|-----|-----|-----|---|
|Set 1|1|2|...|
|Set 2|3|4|...|

To use a lookup table:

* Create a CSV file with your table and add it to your project
* Add the `#meta.inputs`-line to your model and specify a CSV file with a lookup table
* To improve usability, define `caption`, `category` and a tooltip:

```python
#meta.inputs: table {choices: OpenFile("System:AppData/DiffStudio/inputs.csv"); caption: Mode; category: Settings} [Hint]
```

## Performance

Diff Studio solvers ensure fast **in-browser** intergration of ODEs. The following [classic problems](https://archimede.uniba.it/~testset/testsetivpsolvers/?page_id=26#ODE) illustrate their efficiency:

* [Rober](https://archimede.uniba.it/~testset/report/rober.pdf)
  * a stiff system of 3 nonlinear ODEs
  * describes the kinetics of an autocatalytic reaction given by Robertson
* [HIRES](https://archimede.uniba.it/~testset/report/hires.pdf)
  * a stiff system of 8 non-linear equations
  * explains the `High Irradiance Responses' (HIRES) of photomorphogenesis on the basis of phytochrome, by means of a chemical reaction involving eight reactants
* [VDPOL](https://archimede.uniba.it/~testset/report/vdpol.pdf)
  * a system of 2 ODEs proposed by B. van der Pol
  * describes the behaviour of nonlinear vacuum tube circuits
* [OREGO](https://archimede.uniba.it/~testset/report/orego.pdf)
  * a stiff system of 3 non-linear equations
  * simulates Belousov-Zhabotinskii reaction
* [E5](https://archimede.uniba.it/~testset/report/e5.pdf)
  * a stiff system of 4 non-linear ODEs
  * represents a chemical pyrolysis model
* [Pollution](https://archimede.uniba.it/~testset/report/pollu.pdf)
  * a stiff system of 20 non-linear equations
  * describes a chemical reaction part of the air pollution model designed at The Dutch National Institute of Public Health and Environmental Protection

The LSODA, CVODE, MRT, ROS3PRw and ROS34PRw methods demonstrate the following time performance (AMD Ryzen 5 5600H 3.30 GHz CPU):

|Problem|Segment|Points|Tolerance|LSODA, ms|CVODE, ms|MRT, ms|ROS3PRw, ms|ROS34PRw, ms|
|-|-|-|-|-|-|-|-|-|
|[Rober](https://archimede.uniba.it/~testset/report/rober.pdf)|[0, 10E+11]|40K|1E-7|67|85|175|446|285|
|[HIRES](https://archimede.uniba.it/~testset/report/hires.pdf)|[0, 321.8122]|32K|1E-10|125|142|122|362|215|
|[VDPOL](https://archimede.uniba.it/~testset/report/vdpol.pdf)|[0, 2000]|20K|1E-12|268|352|492|1576|760|
|[OREGO](https://archimede.uniba.it/~testset/report/orego.pdf)|[0, 360]|36K|1E-8|76|98|205|483|199|
|[E5](https://archimede.uniba.it/~testset/report/e5.pdf)|[0, 10E+13]|40K|1E-6|7|\*|6|17|8|
|[Pollution](https://archimede.uniba.it/~testset/report/pollu.pdf)|[0, 60]|30K|1E-6|12|15|18|50|23|

\* E5 is skipped for CVODE: the extremely stiff rate constants (spanning 20 orders of magnitude) cause convergence failures with the dense direct linear solver.

This table compares the efficiency of the methods when solving each test problem on a fixed segment and providing solutions at a specified number of points with a given tolerance.

## Platform applications

Once you are satisfied with the result, click <i class="fas fa-layer-plus"></i> **Save to Library** on the top panel to publish your model to the **Library**. It becomes available on the Hub and in **Browse > Apps > Compute > Diff Studio > Library** for yourself and anyone with access to your files.

You can export your model to JavaScript script. To do so:

1. Turn on the **Edit** toggle on the top panel
2. Click `</>` icon. Script editor opens in a new view
3. Click **SAVE** button
4. Script is created, and can be found in the "Scripts" section of the platform

Apply [scripting](https://datagrok.ai/help/compute/scripting) tools to get:

* non-elementary and special functions' use
* Datagrok packages' functions call

## Export to LaTeX and Markdown

Diff Studio can render a model as LaTeX or Markdown for use in papers, reports and notebooks. Click the <i class="fas fa-file-export"></i> **Export** icon on the top panel to open the export dialog. In the dialog:

* pick the output **Format**:
  * `latex` — generates a `.tex` source. Toggle **Standalone** to wrap the output with `\documentclass{article}` and the required `\usepackage` lines so the file compiles with `pdflatex` out of the box
  * `markdown` — generates a `.md` file with LaTeX math blocks, suitable for GitHub, Jupyter, or static-site generators
* review the result in a live preview with syntax highlighting that follows the selected format
* toggle which sections to include (title & description, initial conditions, parameters, constants)
* enable **Compact** mode for small models (no section headings, inline initial conditions)
* choose between `\cdot` and juxtaposition for multiplication
* **Copy** the source to the clipboard or **Download** it as a file

## Citation

If you use Diff Studio in your research, please cite our [JOSS paper](https://doi.org/10.21105/joss.09090):

```
@article{diffstudio2026joss,
  title = {Diff Studio: Ecosystem for Interactive Modeling by Ordinary Differential Equations},
  author = {Viktor Makarichev, Larisa Bankurova, Gennadii Zakharov, Leonid Stolbov, Steven Mehrman, Dan Skatov, Jeffrey Cohen, Paul Sass, Davit Rizhinashvili, Andrew Skalkin},
  journal = {Journal of Open Source Software},
  year = {2026},
  volume = {11},
  number = {120},
  doi = {10.21105/joss.09090}
}
```

## Links

Run Diff Studio online [here](https://public.datagrok.ai/apps/DiffStudio), or complete an interactive [tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations).
Find more features in Diff Studio [docs](https://datagrok.ai/help/compute/diff-studio).

See also

* [Sensitivity analysis](https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis)
* [Parameter optimization](https://datagrok.ai/help/compute/function-analysis#parameter-optimization)
* [Community](https://community.datagrok.ai/t/solving-differential-equations/878)

