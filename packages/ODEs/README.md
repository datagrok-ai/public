# ODEs

ODEs is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform. It provides in-browser tools for solving [intial value problem](https://en.wikipedia.org/wiki/Initial_value_problem) for systems of ordinary differential equations (ODE).

* Go to **Apps** and run **EquaSleek X**
* Enter formulas or modify template
* Press <i class="fas fa-play"></i> **Run** button on the top panel

The solver requires ODEs specification in the declarative form.

## Quick start

Datagrok ensures intuitive tools for the rapid solving ODEs.

* Run the solver: 
  * Go to **Apps** and run **EquaSleek X**
  * Press <i class="fas fa-play"></i> **Run** button on the top panel
* Modify a template:
  * Edit formulas or add new ones
  * Click **F5** or <i class="fas fa-play"></i> **Run** button
* Use an advanced template:
  * Rigth click and select **Templates > Advanced...**
  * Modify formulas and click **F5**
* Save formulas in a local file:
  * Rigth click and select **Save...**
  * Find the file in Downloads
* Load equations from a local file:
  * Rigth click and select **Load...**
  * Select a file with formulas

## Creating project

The platform provides the creation of an ODEs project in a few clicks.

### Templates

Apply built-in templates as a backbone of your project. Rigth click, select **Templates** and choose one the following:

| Template | Features|
|----------|---------|
| `Basic`    | the minimum required |
| `Advanced` | extra math features, including *expressions*, *constants*, *parameters* and *tolerance* specification|
| `Extended` | the *annotating* feature for extended UI generation                 |

### Use cases

The solver has built-in use cases. Get access to them via the context menu. You can use them as a template.

* `Chem reactions`
  * simulates [mass-action kinetics](https://en.wikipedia.org/wiki/Law_of_mass_action)
  * illustrates annotation of inputs
* `Robertson's model`
  * Robertson’s chemical reaction model
  * [stiff equations](https://en.wikipedia.org/wiki/Stiff_equation) example
  * shows how Datagrok solves complicated ODEs
* `Fermentation`
  * models the kinetics of the biochemical reactions in [fermentation](https://en.wikipedia.org/wiki/Fermentation)
  * shows the usage of the `runOnOpen` and  `runOnInput` meta-features
* `PK-PD`
  * simulates pharmacokinetics (PK), pharmacodynamics (PD), and their [relationship](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046)
  * illustrates the usage of the `loop` feature for dosing specification
* `Acid production`
  * models gluconic acid [production](https://oatao.univ-toulouse.fr/9919/1/Elqotbi_9919.pdf) by Aspergillus niger
  * shows the usage of the `update` feature for multi-stage simulation
* `Nimotuzumab`
  * models population pharmacokinetic for [nimotuzumab](https://www.mdpi.com/1999-4923/12/12/1147)
  * demonstrates the `output` feature

Datagrok's ODEs suite has tools for solving both [stiff](https://en.wikipedia.org/wiki/Stiff_equation) and non-stiff equations. It provides a [numerical solution](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations).

### Project structure

A project defines [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem). It contains *name*, *differential equations*, *initial values* and *argument* specifications:

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

You can create extensions by defining components in the control blocks.

Find useful the following:

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
|```#meta.runOnOpen```|Provides computations immediately upon model launch after the project is exported to JavaScript script.|
|```#meta.runOnInput```|Updates results immediately upon input changes after the project is exported to JavaScript script.|

### Platform applications

Datagrok provides an export of each solver project to JavaScript script:

* Press **Export** button on the top panel
* Press **SAVE** button

Find the created script in the platform `Scripts`.

Improve usability. Use `#tags: model` to add your model to `Model Catalog`. Provide a description in the `#description`-line:

```python
#name: Bioreaction
#tags: model
#description: Complex bioreaction simulation
```

Apply annotations in your project to get better UI. Append them, when defining *parameters*, *initial values* and *argument*.

Group inputs by specifying their `category`:

```python
#parameters:
  P1 = 1 {category: Parameters}
  P2 = -1 {category: Parameters}
```

Add `units`:

```python
#inits:
  x = 2 {units: C; category: Initial values}
  y = 0 {units: C; category: Initial values}
```

Define the desired `caption`:

```python
#argument: t
  start = 0 {caption: Initial; units: min; category: Time}
  finish = 2 {caption: Final; units: min; category: Time}
  step = 0.01 {caption: Initial; units: min; category: Time}
```

Provide hints in brackets `[ ]`:

```python
  P1 = 1 {category: Parameters} [P1 parameter tooltip]
```

The export feature provides an extension of your project with [scripting](https://datagrok.ai/help/compute/scripting) tools. You can add the desired functionality:

* results visualization
* non-elementary and special functions use
* Datagrok packages functions call

See also

* [Compute](https://datagrok.ai/help/compute)
* [Scripting](https://datagrok.ai/help/compute/scripting)
* [Viewers gallery](https://datagrok.ai/help/visualize/gallery)
