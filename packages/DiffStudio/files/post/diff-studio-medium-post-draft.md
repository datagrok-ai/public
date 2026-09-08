# Write the Math, Get the App: Solving ODEs in the Browser Without Code

*Diff Studio is both an environment for building ODE models and a hub where they live afterwards. Here is how it works, from the first equation to a shared library.*

---

## The gap between the math and the model

Every modeler knows this moment. The equations are done. They fit on half a page. You know exactly what the system should do.

Then the real work starts: pick a solver, write the integration loop, set up the plotting, wrap it in a notebook. And then a colleague asks, "Can I try it with a different dose?" Now you are either running simulations for them by hand or building a UI you never planned to build.

None of that is the science. It is packaging. And packaging is what quietly decides who gets to use a model: the person who wrote the script, or everyone who has a question about the system it describes. Most of the people who need answers from an ODE model (biologists, clinicians, process engineers, students) never open the script. The model stays with the modeler.

This post is about a different workflow: you write the equations almost the way you would write them on paper, and you get an interactive model in the browser right away. Sliders for every parameter, live plots, fitting to experimental data, and a URL you can send to anyone, whether or not they have ever written code. The point is not to save the modeler ten minutes. The point is to let the model leave the modeler's laptop.

The tool is **Diff Studio**, part of the [Datagrok](https://datagrok.ai) platform. Its solver engine, [Diff Grok](https://github.com/datagrok-ai/diff-grok), is open source (MIT).

*Disclosure: I work at Datagrok. Everything below runs in the free public instance, so you can follow along: [public.datagrok.ai/apps/DiffStudio](https://public.datagrok.ai/apps/DiffStudio).*

![run](./run.gif)

Here is what we will build up to:

- a working model from a dozen lines of plain text that read like math,
- parameters that automatically become UI controls,
- a two-compartment pharmacokinetic model,
- fitting, sensitivity analysis, and sharing: the parts that are painful to bolt onto a script,
- how models accumulate into a library your whole team can browse and reuse.

Let's start with the smallest possible example.

---

## Your first model in 60 seconds

Open Diff Studio (**Apps > Diff Studio**), pick the **Basic** template or just switch on the **Edit** toggle, and replace the text with this classic predator–prey system:

```
#name: Lotka-Volterra

#equations:
  dx/dt = alpha * x - beta * x * y
  dy/dt = delta * x * y - gamma * y

#inits:
  x = 10
  y = 5

#argument: t
  initial = 0
  final = 50
  step = 0.05

#parameters:
  alpha = 1.1
  beta = 0.4
  delta = 0.1
  gamma = 0.4
```

Press **F5** (or click **Refresh**). That's it. The model is solved and plotted.

![create](./create.gif)

A quick tour of what each block does:

- **`#name`** — the model identifier. It shows up in the tab and, later, in the model library.
- **`#equations`** — the system itself. Write derivatives as `dx/dt`, use any variable names you like (multi-letter names are fine), and use ordinary math notation: `*`, `/`, `exp()`, `cos()`, and so on.
- **`#inits`** — initial conditions for every function you are solving.
- **`#argument`** — the independent variable, its range, and the output step. The solver returns values on this grid.
- **`#parameters`** — quantities that appear in the equations and that you want to play with. Each one becomes an input in the UI. (Values that should stay fixed go under `#constants`; more on that below.)

There is nothing else to configure. No solver selection, no tolerances, no plotting code. The default solver handles both stiff and non-stiff systems; you can override it later if you need to.

Two small things that matter more than they seem. You can add comments with `//` right next to an equation. And the whole model is plain text: download it as an `.ivp` file, edit it in any editor, keep it in git, paste it into a chat.

---

## The same model in Python, for an honest comparison

To be clear about what Diff Studio is and isn't replacing, here is the equivalent in Python with SciPy, the tool most of us would reach for:

```python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

alpha, beta, delta, gamma = 1.1, 0.4, 0.1, 0.4

def rhs(t, z):
    x, y = z
    return [alpha * x - beta * x * y,
            delta * x * y - gamma * y]

t_eval = np.linspace(0, 50, 1001)
sol = solve_ivp(rhs, (0, 50), [10, 5], t_eval=t_eval)

plt.plot(sol.t, sol.y[0], label="x (prey)")
plt.plot(sol.t, sol.y[1], label="y (predator)")
plt.xlabel("t"); plt.legend(); plt.show()
```

The Python version is short, readable, and perfectly good. SciPy is excellent, and if you need a custom pipeline around your solver, Python (or Julia) is the right place to be.

But look at what the script is made of. Two lines are the actual model. The rest is glue: unpacking the state vector, building the time grid, calling the integrator, labeling the plot. And the script still does nothing the Diff Studio version doesn't already do.

Now answer the colleague's question, "what if `beta` is 0.6?", with each version. In Python you edit a number and rerun. In Diff Studio you drag a slider. That difference is the whole point, and the next section is about making the sliders good.

---

## Parameters become the interface

The model already has four inputs, but they are bare number fields labeled with variable names. Diff Studio builds the UI from annotations you attach to any input, in curly braces after the value. Let's make the Lotka–Volterra model something a non-modeler could pick up:

```
#parameters:
  alpha = 1.1 {min: 0; max: 3; step: 0.05; caption: Prey growth rate; category: Prey}
  beta  = 0.4 {min: 0; max: 2; step: 0.05; caption: Predation rate;   category: Prey}
  delta = 0.1 {min: 0; max: 1; step: 0.01; caption: Predator growth;  category: Predator}
  gamma = 0.4 {min: 0; max: 2; step: 0.05; caption: Predator death;   category: Predator}

#inits:
  x = 10 {min: 1; max: 50; caption: Initial prey;     category: Population; units: ind}
  y = 5  {min: 1; max: 50; caption: Initial predator; category: Population; units: ind} [Number of predators at t = 0]
```

Click the Edit toggle.

![annotate](./annots.gif)

What each annotation does:

- **`min` / `max`** — turns the field into a slider. Add **`step`** to control the increment.
- **`caption`** — the label people see instead of `alpha`.
- **`category`** — groups related inputs into collapsible sections.
- **`units`** — shown next to the field. Cosmetic, but it prevents a lot of "is this in hours or minutes?" questions.
- **`[ ... ]`** — a tooltip. Anything in square brackets after the annotation appears on hover.

The same annotations work in the `#argument` block, so `initial`, `final`, and `step` can get captions and ranges too.

Two more things worth knowing here:

**`#parameters` vs `#constants`.** Both are named values you can use in equations. The difference is only what shows up in the UI. If a value is a physical constant or something users should never touch, move it to `#constants` and the panel stays uncluttered:

```
#constants:
  K = 100   // carrying capacity, fixed for this study
```

**Sliders are live.** Every time a slider moves, the system is re-solved and the plot updates. For a model like this one the update is effectively instant, which is what makes exploration feel different from "edit, rerun, look".

You now have an interactive predator–prey app that took a couple of minutes and contains no code, only equations and labels. Anyone who understands what "prey growth rate" means can use it, and that is a much larger group than the one that can edit a `solve_ivp` call. Let's move to something closer to real work.

---

## A real example: two-compartment pharmacokinetics

Diff Studio ships with a **Library** of ready models: PK and PK-PD simulations, a bioreactor, fermentation, chemical kinetics, and a few classic stiff benchmarks. The model below is a compact variant of the pharmacokinetics example there; you can paste it in as is.

A two-compartment model with first-order absorption from the gut:

```
#name: Two-compartment PK

#equations:
  dA_gut/dt  = -ka * A_gut
  dA_c/dt    =  ka * A_gut - (CL / Vc) * A_c - Q * (A_c / Vc - A_p / Vp)
  dA_p/dt    =  Q * (A_c / Vc - A_p / Vp)

#expressions:
  C_central    = A_c / Vc
  C_peripheral = A_p / Vp

#inits:
  A_gut = 500 {caption: Dose; units: mg; category: Dosing; min: 0; max: 1000}
  A_c   = 0   {caption: Central amount;    units: mg; category: Initial state}
  A_p   = 0   {caption: Peripheral amount; units: mg; category: Initial state}

#parameters:
  ka = 1.0 {caption: Absorption rate;     units: 1/h;  category: PK parameters; min: 0.1; max: 5}
  CL = 5.0 {caption: Clearance;           units: L/h;  category: PK parameters; min: 0.5; max: 20}
  Vc = 30  {caption: Central volume;      units: L;    category: PK parameters; min: 5;   max: 100}
  Vp = 60  {caption: Peripheral volume;   units: L;    category: PK parameters; min: 5;   max: 200}
  Q  = 3.0 {caption: Intercompartmental clearance; units: L/h; category: PK parameters; min: 0.1; max: 20}

#argument: t
  initial = 0  {caption: Start; units: h}
  final   = 48 {caption: End;   units: h}
  step    = 0.1 {caption: Step; units: h}

#output:
  t            {caption: Time, h}
  C_central    {caption: Central conc., mg/L}
  C_peripheral {caption: Peripheral conc., mg/L}
```

You get:

![pk](./pk.png)

Two new blocks appear here.

**`#expressions`** are auxiliary quantities computed from the state, parameters, and time, with no ODEs involved. Pharmacokinetics is the perfect illustration: the equations are naturally written in *amounts*, but nobody plots amounts. They plot *concentrations*. Instead of rewriting the system, you define `C_central = A_c / Vc` once and use it anywhere, including in the output.

**`#output`** controls what ends up in the results table and on the chart. By default, Diff Studio shows every solved function. Here we hide the raw amounts and show only what a pharmacologist wants to see, with proper column headers. Anything from `#equations` or `#expressions` can go here.

Now the colleague's question, "what happens with a 250 mg dose in a patient with half the clearance?", is two slider moves. And they can do it themselves.

Everything so far could, with some effort, be reproduced in a notebook. The next section is about the parts that usually can't, at least not without building a small product around your script.

---

## The parts that are hard to bolt onto a script

A model in a notebook answers the questions its author thought of. A model in Diff Studio comes with three tools that answer the questions everyone else asks.

### Fit: "Which parameters explain my data?"

You have measured concentrations at a few time points and want the PK parameters that reproduce them. Click the **Fit** icon on the top panel, load your table with the observed values, pick which parameters are allowed to vary and in what ranges, and run.

Diff Studio searches the parameter space, reports the goodness of fit, and overlays the fitted curve on the data. The optimization runs in parallel in the browser, so even complex models come back quickly.

![fit](./fit.gif)

Nothing changed in the model text. Fitting is a property of *any* Diff Studio model, not something you implement per model.

### Sensitivity: "Which parameters actually matter?"

Click **Sensitivity**. Choose the inputs to vary and the outputs to watch, pick a method (Monte Carlo, Sobol, or Grid), and run. You get a picture of how each parameter drives each output, which is usually the first thing a reviewer asks about a model.

![sa](./sa.gif)

Again: no extra code. This is the same panel for the predator–prey toy and for a 20-equation bioreactor.

### Share: "Can I try it?"

A model run, including its parameter values, is encoded in the URL. Copy the link from the address bar and send it. The recipient opens the exact run you were looking at, with the same sliders, and can start changing things immediately. No environment to install, no notebook to re-execute.

When you need the model in a paper or a report, the **Download** menu exports it to Markdown or LaTeX with the equations properly typeset.

Together, these three cover most of "can you rerun it with…", "how sensitive is it to…", and "can I have a copy…": the requests that quietly turn a two-week modeling task into a two-month one.

---

## From one model to a hub

Sharing a link solves the problem for one model and one colleague. The larger problem is that in many organizations ODE models live in scattered scripts, notebooks, and spreadsheets, each with its own author, its own conventions, and its own way of running. Nobody knows which model of the bioreactor is current, and the PK model from last year's project has to be rediscovered from scratch.

Diff Studio is built to be the place where those models live. Click **Save to Library** and the model becomes part of a catalog that is searchable and shared under the platform's access rules. Attach a help page with the assumptions and references. Anyone with access opens it from the same **Library** tab we used earlier, runs it, fits it to their data, and never needs to ask the author for the file.

![hub](./hub.png)

A model doesn't have to stay inside Diff Studio, either. Convert it to a Datagrok script with one click and it becomes a function like any other on the platform, a building block for pipelines, dashboards, and custom applications: a dosing calculator built on the PK model, a what-if tool for the process engineers, a teaching demo. The equations stay the same; what changes is who can reach them.

This is why we describe Diff Studio as an environment *and* a hub: an environment for building and exploring ODE models without code, and a hub where a team's models are collected, shared, and reused instead of rewritten.

---

## What's under the hood

Diff Studio's engine is [Diff Grok](https://github.com/datagrok-ai/diff-grok), an open-source TypeScript library (MIT license, zero dependencies) for initial value problems. It solves both stiff and non-stiff systems directly in the browser, without a server round-trip.

The default solver is a Rosenbrock–Wanner method (ROS34PRw), which handles stiff problems like the Robertson and HIRES benchmarks without the user having to know they are stiff. If you want control, the `#meta.solver` block lets you switch to other methods (explicit Runge–Kutta schemes, Adams multistep methods, or the adaptive LSODA and CVODE) and set time limits and tolerances:

```
#meta.solver: {method: 'lsoda'; maxTimeMs: 100}
#tolerance: 0.00001
```

For readers who want the details (the numerical methods, the computational pipeline that makes in-browser solving fast, and benchmark results against reference solvers), the design is described in two papers:

- **Journal of Open Source Software (2026):** [Diff Studio: Ecosystem for Interactive Modeling by Ordinary Differential Equations](https://doi.org/10.21105/joss.09090) — covers Diff Grok, the open-source engine.
- **Springer, CoMeSySo 2025 proceedings:** [Diff Studio: Web-Based Environment for Interactive Modeling with Ordinary Differential Equations](https://doi.org/10.1007/978-3-032-22236-7_31) — the foundational paper on the web-based approach and its performance on classic benchmark problems.

Stiff systems in the browser deserve a post of their own. That one is coming.

---

## Try it

The whole idea in one sentence: **you write the math, and the platform builds the app around it.** Equations become an interactive model, parameters become sliders, and fitting, sensitivity analysis, and sharing come for free. Models stop being personal scripts and become shared assets that anyone on the team can open, question, and build on.

Scientific computing has spent decades getting more powerful and, at the same time, more concentrated in the hands of people who write code. We think the next step goes the other way: putting well-built models within reach of the people who have the questions. Diff Studio is one attempt at that.

Three ways to start:

- **Run it now:** [Diff Studio on public.datagrok.ai](https://public.datagrok.ai/apps/DiffStudio). Free, nothing to install. Start from a template or open a model from the Library.
- **Take the guided tour:** the [interactive tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations) inside the platform walks you through building and exploring a model step by step.
- **Read the docs:** [Diff Studio documentation](https://datagrok.ai/help/compute/diff-studio) has the full syntax reference, and the [community thread](https://community.datagrok.ai/t/solving-differential-equations/878) tracks new features as they land.

If there is a specific model you'd like to see built this way (a SIR epidemic, enzyme kinetics, a PK-PD system with an effect compartment), say so in the comments. The next posts in this series will follow what readers ask for.
