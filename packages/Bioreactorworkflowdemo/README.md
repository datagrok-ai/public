# Bioreactor Workflow Demo

`Bioreactorworkflowdemo` is a Datagrok package that demonstrates the main
low-code workflow capabilities of Compute2 through a simple bacterial growth
simulation. Users configure a bioreactor, add and edit daily cultivation
steps, propagate biomass and substrate between days, validate the experiment,
and collect the completed run in a summary.

## Workflow overview

```text
Bioreactor configuration
        | initial biomass, substrate, and model parameters
        v
Daily cultivation (dynamic workflow)
        |-- Day 1
        |-- Day 2       add, remove, or reorder days
        `-- Day N
        |
        v
Review and summary
```

The daily calculations are placed in a nested dynamic workflow between two
fixed steps. Each day consumes the final state of the preceding day and emits
the starting state for the next one. Day 1 receives its initial state from the
bioreactor configuration. The final fixed step aggregates all daily outputs
into the experiment summary.

## Compute2 capabilities demonstrated

| Workflow component | Compute2 capability |
| --- | --- |
| Bioreactor configuration | Static script step with an automatically generated input form |
| Daily cultivation | Nested dynamic workflow using `stepTypes` and `initialSteps` |
| Day-to-day state propagation | Reactive data links and consistency restrictions |
| Input limits | IO validators with errors, warnings, and suggested corrections |
| Number and order of days | A pipeline validator that reacts to structural changes |
| Conditional daily fields | Metadata links for choices and conditionally visible inputs |
| Dynamic day names | Node metadata such as `Day 3 - Biomass 4.2 g/L` |
| User-triggered operations | Actions for correcting inputs or adding predefined stages |
| Final results | A summary step that receives outputs from all dynamic days |
| Saving results | Compute2 history and a custom Excel or CSV export |

## Bioreactor model

The initial demonstration uses a deterministic, discrete daily model. It tracks:

- biomass concentration (`X`, g/L);
- substrate concentration (`S`, g/L);
- reactor volume;
- growth rate based on Monod kinetics;
- biomass decay;
- substrate consumption based on biomass yield;
- optional feed and harvest operations.

A daily step accepts cultivation conditions and the incoming reactor state,
then calculates the final biomass, substrate, and volume. Reordering the steps
therefore represents changing the order of cultivation periods rather than
changing immutable calendar dates.

## Validation

The workflow is intended to demonstrate both input and structural validation:

- initial biomass, reactor volume, cultivation duration, and kinetic constants
  must be positive;
- carrying capacity must exceed the initial biomass;
- feed and harvest operations must leave a valid reactor volume;
- calculated biomass and substrate concentrations cannot be negative;
- the workflow must contain at least one daily cultivation step;
- the number of daily steps is limited to a configurable demonstration maximum;
- unusual pH or temperature and near-depleted substrate produce warnings.

For clarity and responsive interaction, the demo is intended for approximately
one to fourteen daily steps. The summary calculation remains pure, while saving
or exporting the completed experiment is implemented as a separate action.
