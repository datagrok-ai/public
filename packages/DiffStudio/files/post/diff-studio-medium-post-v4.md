# Write the math, get the app: solving ODEs in the browser without code

*By Viktor Makarichev, Scientific Application Developer at Datagrok, Inc.*

---

The math behind a two-compartment PK model is three equations and five parameters. Exposing it to users, however, is a different job entirely.

Someone has to choose and configure an ODE solver. Someone has to build the plots. Someone has to turn parameters into inputs, decide which ones users can change, add units and labels, handle repeated dosing, save and load runs, and connect the model to experimental data. Then someone has to maintain all of that code when the model changes.

The model, in other words, is often the small part of building a usable modeling application. The scaffolding around it is the work.

That work is unavoidable. ODE-based models are used throughout biopharma R&D: pharmacokinetics and pharmacodynamics, quantitative systems pharmacology, bioprocess modeling, and more. These models let scientists explore systems that are expensive or impossible to experiment on directly. A simulation is not a substitute for an experiment, but it can make the space of plausible decisions much easier to explore.

So what options are there?

The usual way to put a model in front of people is to write it in R, Python, MATLAB, or another language, then build a notebook, dashboard, or web application around it. That is a perfectly reasonable approach, especially when the application needs highly customized behavior.

But there is a recurring problem.

A parameter that was not exposed requires a code change. A different output requires a code change. A new UI control requires a code change. A modification to the model means updating the implementation around it.

But what if the model could describe its own interface?

---

## The idea behind Diff Studio

That is the idea behind **Diff Studio**, part of the [Datagrok](https://datagrok.ai) platform.

Instead of writing the model as application code, you describe the math directly:

```
#equations:
  dA_gut/dt = -ka * A_gut
  dA_c/dt   = ka * A_gut - (CL / Vc) * A_c
```

Then define the initial conditions and parameters:

```
#inits:
  A_gut = 500 {caption: Dose; units: mg; min: 0; max: 1000}

#parameters:
  CL = 5 {min: 0.5; max: 20; caption: Clearance; units: L/h}
```

Those annotations are not application code. They are part of the model description, but they give the platform enough information to generate the corresponding interface: a labeled input, a range, a unit, a slider.

The important part is not that a slider appears.

The important part is that nobody wrote the slider.

Diff Studio solves the equations, generates the interface, and updates the visualization when an input changes. The model remains a mathematical description rather than becoming another piece of custom software.

You write the equations. The platform handles the machinery around them.

---

## A PK model is a good test

Consider a two-compartment pharmacokinetic model with first-order absorption (the dose is assumed to be fully absorbed):

```
#name: Two-compartment PK

#equations:
  dA_gut/dt = -ka * A_gut
  dA_c/dt   = ka * A_gut
              - (CL / Vc) * A_c
              - Q * (A_c / Vc - A_p / Vp)
  dA_p/dt   = Q * (A_c / Vc - A_p / Vp)

#expressions:
  C_central    = A_c / Vc
  C_peripheral = A_p / Vp

#inits:
  A_gut = 500 {caption: Dose; units: mg; min: 0; max: 1000}
  A_c   = 0
  A_p   = 0

#parameters:
  ka = 1.0 {caption: Absorption rate; units: 1/h; min: 0.1; max: 5}
  CL = 5.0 {caption: Clearance; units: L/h; min: 0.5; max: 20}
  Vc = 30  {caption: Central volume; units: L; min: 5; max: 100}
  Vp = 60  {caption: Peripheral volume; units: L; min: 5; max: 200}
  Q  = 3.0 {caption: Intercompartmental clearance; units: L/h; min: 0.1; max: 20}

#argument: t
  initial = 0
  final = 48
  step = 0.1
```

There is a lot of information here, but very little "software code". The equations describe the system; the parameters describe what can vary; the annotations describe how those inputs should appear to a user.

**`#expressions`** lets you define quantities derived from the solved state—in this case, concentrations from compartment amounts—and **`#output`** can control what appears in the results.

Change clearance and the model is solved again. Change the dose and it is solved again. Change the volume or the absorption rate and the plots update.

You are exploring the model rather than editing an application.

---

## The application layer becomes a platform capability

Once the model is represented separately from the surrounding application, functionality that would normally be implemented repeatedly can operate on the model itself.

Take fitting.

Suppose you have observed concentration data and want to find parameter values that reproduce them. In a conventional application, fitting can become another piece of model-specific plumbing: define the objective, connect it to the solver, select the parameters, constrain the search, and display the result.

In Diff Studio, fitting is an operation on the model.

Sensitivity analysis works the same way. You can vary inputs and examine their influence on outputs using methods including Monte Carlo sampling and Sobol sensitivity indices.

The reusable object is not the application. It is the model.

Once you have the model, simulation, fitting, sensitivity analysis, visualization, and sharing become just some of the things you can do with it.

---

## The model can also become a shared artifact

There is another problem with scientific models that has nothing to do with differential equations.

Models live in notebooks, scripts, project folders, shared drives, and elsewhere. Duplicates or slightly different versions appear across projects, and finding or maintaining the right model can become difficult.

Datagrok treats the model as something that can live in a shared environment.

A model can be saved to the library, where it can be shared and reopened by other users. Individual model runs can also be shared by URL. The model itself remains editable text and can be downloaded or exported to formats such as Markdown and LaTeX.

That matters because collaboration around models is often just as important as running them.

A scientist may want to send a colleague a particular parameterization, not merely the equation.

A model developer may want to publish the assumptions alongside the model.

A team may want one place where its PK, PK-PD, bioreactor, and kinetic models live instead of maintaining separate collections of scripts.

The interface for those activities does not have to be reinvented for every model.

---

## From model to application

There is a useful middle ground between "no code" and "build everything from scratch."

How much custom software should you have to build around mathematics before another scientist can use the model?

For many models, the answer does not need to be very much.

A declarative model lets the platform take responsibility for the repetitive parts: interpreting the equations, generating inputs, solving the system, plotting the results, and exposing common analysis operations.

Under the hood, Diff Studio uses [Diff Grok](https://github.com/datagrok-ai/diff-grok), an open-source engine for solving initial value problems for ordinary differential equations. It supports both stiff and non-stiff systems and multiple numerical methods — Rosenbrock–Wanner, Runge–Kutta and Adams families — as well as the LSODA and CVODE solvers with automatic stiffness detection.

There is also a path from an interactive model to a custom application. A Diff Studio model can be converted into a script, preserving its input annotations. The same model can then move into a larger workflow or a specialized application when more customization is needed.

So the progression does not have to be:

```
equation → rewrite as software → build application
```

It can be:

```
equation → interactive model → shared model → specialized application
```

A researcher starts by asking, "What happens if clearance is lower?"

Then, "Can we fit these observations?"

Then, "Can everyone on the project use the same model?"

Then, "Can we put it into our workflow?"

Those are different questions, but they all start with the same mathematical object.

---

## Try it

Developing a good model is hard. The science can take years.

But the mathematical and numerical methods for solving ODEs have been around for over a century. What gets rebuilt for every new model is the software around it: the interface, plots, fitting, analysis, sharing, and maintenance.

That is the problem Diff Studio solves.

The model remains a model. The platform turns it into something people can use.

**Write the math. Get the app.**
