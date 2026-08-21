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

The initial demonstration uses a well-mixed fed-batch model with one
growth-limiting substrate. The user adds material once, at the beginning of
each day, and the culture then grows for 24 hours without further inflow or
outflow.

### Assumptions

- `solutionAdded` is the volume of a substrate-free solution, such as water or
  buffer. It contains no biomass and all nutrients other than the modeled
  substrate are assumed non-limiting.
- `substrateAdded` is the dry mass of the feeding substrate. Its volume is
  negligible and it dissolves instantaneously.
- The reactor is ideally and instantaneously mixed after each addition.
- Temperature, pH, and oxygen are constant and non-limiting, so their effects
  are included implicitly in the kinetic constants.
- Biomass decay does not return reusable substrate to the medium.
- Sampling, evaporation, product formation, inhibition, and harvest are not
  included in this version.

These definitions are important. If the substrate is supplied as a liquid feed
rather than as a dry mass, use `feedVolume` and `feedSubstrateConcentration` as
the daily inputs instead. In that formulation,
`substrateAdded = feedVolume * feedSubstrateConcentration`; adding a separate
solution volume as well would represent a second, substrate-free stream.

### Symbols and configuration parameters

| Symbol / parameter | Type | Units | Valid values | Description |
| --- | --- | --- | --- | --- |
| `V0` | double | L | `V0 > 0` | Initial liquid volume |
| `X0` | double | g biomass/L | `X0 >= 0` | Initial biomass concentration |
| `S0` | double | g substrate/L | `S0 >= 0` | Initial limiting-substrate concentration |
| `muMax` (`mu_max`) | double | 1/h | `muMax > 0` | Maximum specific growth rate |
| `monodKs` (`K_s`) | double | g substrate/L | `monodKs > 0` | Substrate concentration at which the Monod growth rate is half of `muMax` |
| `biomassYield` (`Y_XS`) | double | g biomass/g substrate | `biomassYield > 0` | Biomass formed per unit of consumed substrate |
| `decayRate` (`k_d`) | double | 1/h | `decayRate >= 0` | First-order biomass decay rate |
| `dayDuration` (`Delta t`) | fixed double | h | `24` | Duration simulated by every daily step |
| `maxVolume` | double | L | `maxVolume > V0` | Reactor working-volume limit used for validation |

The configuration step produces the initial reactor state
`(V_0, X_0, S_0)`. Kinetic parameters remain fixed for all days in one
workflow run.

### Per-day inputs

Only the first two parameters are editable by the user. The remaining values
are read-only inputs propagated from the configuration or preceding day.

| Parameter | Type | Units | Valid values | Description |
| --- | --- | --- | --- | --- |
| `solutionAdded` (`Delta V_n`) | double | L | `Delta V_n >= 0` | Substrate-free solution added at the start of day `n` |
| `substrateAdded` (`Delta M_S,n`) | double | g substrate | `Delta M_S,n >= 0` | Limiting substrate added at the start of day `n` |
| `incomingVolume` (`V_n-1`) | linked double | L | `V_n-1 > 0` | Final liquid volume from the preceding day |
| `incomingBiomass` (`X_n-1`) | linked double | g biomass/L | `X_n-1 >= 0` | Final biomass concentration from the preceding day |
| `incomingSubstrate` (`S_n-1`) | linked double | g substrate/L | `S_n-1 >= 0` | Final substrate concentration from the preceding day |

For Day 1, the incoming values are `V0`, `X0`, and `S0` from the
configuration step.

### Addition and dilution formulas

For day `n`, convert the incoming concentrations to masses and apply the two
additions:

```text
M_X,in = X_n-1 * V_n-1
M_S,in = S_n-1 * V_n-1

V_n(0) = V_n-1 + Delta V_n
M_X,n(0) = M_X,in
M_S,n(0) = M_S,in + Delta M_S,n

X_n(0) = M_X,n(0) / V_n(0)
S_n(0) = M_S,n(0) / V_n(0)
```

Adding substrate-free solution therefore dilutes both biomass and residual
substrate concentrations but does not change either mass. Adding dry substrate
increases substrate mass without changing the volume.

### Growth formulas

During the following 24-hour cultivation period, volume is constant and the
specific growth rate follows Monod kinetics:

```text
mu(S) = mu_max * S / (K_s + S), for S > 0
mu(0) = 0

dX/dt = (mu(S) - k_d) * X
dS/dt = -mu(S) * X / Y_XS
dV/dt = 0
```

The equations are solved on `0 <= t <= 24 h` from the post-addition state
`(X_n(0), S_n(0), V_n(0))`. Substrate uptake stops at `S = 0`; from that point
`mu = 0`, `S` remains zero, and biomass changes only through decay. Numerical
roundoff must never be allowed to produce negative biomass or substrate.

This coupled Monod system is integrated with adaptive RK4 step-doubling using
a relative and absolute tolerance of `1e-8`; treating the start-of-day growth
rate as constant for the whole day would overestimate growth when substrate
becomes depleted. The daily profile contains hourly samples and an additional
point at substrate depletion when depletion occurs between samples.

### Per-day outputs

| Parameter | Type | Units | Formula / description |
| --- | --- | --- | --- |
| `finalVolume` (`V_n`) | double | L | `V_n(24)`; equals `V_n-1 + Delta V_n` |
| `finalBiomass` (`X_n`) | double | g biomass/L | `X_n(24)` |
| `finalSubstrate` (`S_n`) | double | g substrate/L | `S_n(24)` |
| `finalBiomassMass` (`M_X,n`) | double | g biomass | `X_n * V_n` |
| `finalSubstrateMass` (`M_S,n`) | double | g substrate | `S_n * V_n` |
| `biomassMassChange` | double | g biomass | `M_X,n - M_X,in` |
| `substrateConsumed` | double | g substrate | `M_S,in + Delta M_S,n - M_S,n` |
| `dailyProfile` | dataframe | mixed | Time, volume, biomass concentration, substrate concentration, and their masses during the day |

`finalVolume`, `finalBiomass`, and `finalSubstrate` become the linked incoming
state of the next dynamic step. The summary step combines every day's inputs,
outputs, and `dailyProfile` into the complete cultivation history. Reordering
daily steps represents changing the order of cultivation periods rather than
changing immutable calendar dates.

The model follows the conventional Monod relationship and fed-batch biomass,
substrate, and volume balances used in microbial cultivation models. See the
[experimental evaluation of Monod kinetics](https://doi.org/10.1128/am.22.6.1041-1047.1971),
a [standard fed-batch balance example](https://pmc.ncbi.nlm.nih.gov/articles/PMC11630647/),
and a [pulse-fed cultivation model](https://pmc.ncbi.nlm.nih.gov/articles/PMC6999281/).

## Validation

The workflow is intended to demonstrate both input and structural validation:

- initial volume, `muMax`, `monodKs`, `biomassYield`, and `maxVolume` must
  be positive;
- initial biomass, initial substrate, `decayRate`, and daily additions cannot
  be negative;
- the volume after a daily addition cannot exceed `maxVolume`;
- calculated biomass and substrate concentrations cannot be negative;
- the workflow must contain at least one daily cultivation step;
- the number of daily steps is limited to a configurable demonstration maximum;
- substrate depletion before the end of a day, excessive dilution, or
  `decayRate >= muMax` produces a warning.

For clarity and responsive interaction, the demo is intended for approximately
one to fourteen daily steps. The summary calculation remains pure, while saving
or exporting the completed experiment is implemented as a separate action.
