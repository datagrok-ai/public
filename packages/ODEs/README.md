# Diff Studio

Diff Studio is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai) platform. It provides in-browser tools for solving [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem) for systems of ordinary differential equations (ODE) and its interactive exploration.

* Go to **Apps** and run **Diff Studio**
* Enter formulas or modify template
* Click **F5** or go to **Run** tab

The solver takes a set of the differential equations in a declarative form, and creates a UI that solves the equations, visualizes the results, and lets you change parameters on the fly.

## Quick start

Datagrok provides intuitive tools for the rapid solving ODEs.

* Run the solver:
  * Go to **Apps** and run **Diff Studio**
  * Go to **Run** tab to launch computations
* Modify a template:
  * Edit formulas or add new ones
  * Click **F5** or go to **Run** tab
  * Cahnge inputs and explore your model
* Use an advanced template:
  * Press **Open** icon on the top panel
  * Select **Templates > Advanced...**
  * Modify formulas and click **F5**
* Save formulas in a local file:
  * Press **Save** icon on the top panel  
  * Find the *ivp*-file in Downloads, modify it using any text editor
* Drag-n-drop:
  * Drag *ivp*-file with equations right into the browser
* Load equations from a local file:
  * Press **Open** icon on the top panel
  * Select **From file...**
  * Select a file with formulas

## Create model from template

Start from one of these templates:

| Template | Features|
|----------|---------|
| `Basic`    | the simplest model|
| `Advanced` | extra math features, including *expressions*, *constants*, *parameters* and *tolerance* specification|
| `Extended` | the *annotating* feature for extended UI generation                 |

## Use cases

The solver has built-in use cases. Get access to them via the context menu. You can use them as a template.

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
  * illustrates the capability of modeling complex processes described by a large number of equations

Datagrok's ODEs suite has tools for solving both [stiff](https://en.wikipedia.org/wiki/Stiff_equation) and non-stiff equations. It provides a [numerical solution](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations).

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
|```#name```|Defines a name.|
|```#equations```|Differential equation specification|
|```#argument```|Independent variable specification|
|```#inits```|Initial values specification|
|```#constants```|Constants specification|
|```#parameters```|Parameters specification|
|```#expressions```|Defines auxiliary compuations.|
|```#output```|Defines output columns and their captions|
|```#tags```|Specifies tags (`model`, `app`, etc.).|
|```#description```|Defines description of the model.|
|```#comment```|Specifies comments block.|

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

## Export to JavaScript script

Once you are satisfied with the result, you can convert your model to a Datagrok application. To do so:

1. Press **JS** button on the top panel
2. Press **SAVE** button
3. Script is created, and can be found in the "Scripts" section of the platform

Improve usability. Use `#tags: model` to add your model to `Model Catalog`. Provide a description in the `#description`-line:

```python
#name: Bioreaction
#tags: model
#description: Complex bioreaction simulation
```

Apply [scripting](https://datagrok.ai/help/compute/scripting) tools to get:

* non-elementary and special functions' use
* Datagrok packages' functions call

Find more features in Diff Studio [docs](https://datagrok.ai/help/compute/diff-studio).

See also

* [Sensitivity analysis](https://datagrok.ai/help/compute/#sensitivity-analysis)
* [Viewers gallery](https://datagrok.ai/help/visualize/gallery)
* [Scripting](https://datagrok.ai/help/compute/scripting)
