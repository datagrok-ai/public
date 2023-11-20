---
title: "From scratch"
sidebar_position: 2
---

Explore built-in templates of the [solver](solver.md). They cover all capabilities. Use them as a backbone of custom projects.

| Template   | Comment                                   |
|------------|-------------------------------------------|                                
| `Basic`    | the minimum required project              |
| `Advanced` | extra math features                       |
| `Extended` | [extensions](extensions.md) design sample |

Right-click and select **Templates** to get one of them

![add-to-workspace](from-scratch.gif)

## Project structure

A project defines [initial value problem](https://en.wikipedia.org/wiki/Initial_value_problem). It contains *name*, *differential equations*, *initial values* and *argument* specifications.

Start with the line defining a name:

```
#name: Problem1
```

Place differential equations in the `#equations`-block:

```
#equations:
  dx/dt = x + y + exp(t)
  dy/dt = x - y - cos(t)
```

Set variable, its *initial*, *final* and *step* in the `#argument`-block:

```
#argument: t
  initial = 0
  final = 1
  step = 0.01
```

Datagrok provides a numerical solution within the range *[initial, final]* with the specified grid *step*.

Define initial values in the `#inits`-block:

```
#inits:
  x = 2
  y = 5
```

## Advances

The solver supports the use of *constants*, *parameters* and *expressions*.

Specify constants in the `#constants`-block

```
#constants:
  C1 = 1
  C2 = 3
```

Set parameters in the `#parameters`-block:

```
#parameters:
  P1 = 1
  P2 = -1
```

Define auxiliary compuations in the `expressions`-block:

```
#expressions:
  E1 = C1 * t + P1
  E2 = C2 * cos(2 * t) + P2
```

Expressions may depend on constants, parameters and the argument. Equations may contain results of their evaluation.

Find more capabilities in

* [Extensions](extensions.md)
