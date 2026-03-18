# SirSeirEpidemic

`SirSeirEpidemic` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform

This package is created using the following CC skills:

1. `create-interactive-scientific-application-spec`

User input:
```bash
/create-interactive-scientific-application-spec use the prompt from PROMPT.md. Save result to SPEC.md
```

Output: [SPEC.md](SPEC.md)

2. `implement-interactive-scientific-application-from-spec`

Manual updates:

* bumped version of `diff-grok`


## Prompt: SIR / SEIR Epidemic Model

Create spec of an interactive application of the SIR and SEIR epidemiological models.

### Mathematical Model

#### SIR:
```
dS/dt = -β·S·I/N
dI/dt = β·S·I/N - γ·I
dR/dt = γ·I
```

#### SEIR (extended, toggle to switch):
```
dS/dt = -β·S·I/N
dE/dt = β·S·I/N - σ·E
dI/dt = σ·E - γ·I
dR/dt = γ·I
```

Where N = S + I + R (or S + E + I + R for SEIR).

Default parameters: β = 0.3, γ = 0.1 (R₀ = β/γ = 3), σ = 0.2 (for SEIR), N = 10,000.
Initial conditions: S₀ = N − 1, I₀ = 1, R₀ = 0.

### ODE Solver

- Use the **MRT (Modified Rosenbrock Triple)** method from the **diff-grok** library for numerical integration.
- Step size: adaptive, initial h = 0.1.
- Integration interval: t ∈ [0, 300] days.

### Visualization

#### Line Chart — main (epidemic dynamics)
- **Line chart** with curves S(t), E(t) (if SEIR), I(t), R(t).
- Colors: S — blue, E — orange, I — red, R — green.
- Area under I(t) lightly filled red (area chart) to emphasize the infection peak.
- Vertical dashed line at the peak of I(t) with label "Peak: day N, I = X".

#### Line Chart — secondary (effective reproduction number)
- **Line chart**: R_eff(t) = R₀ · S(t)/N.
- Horizontal line at R_eff = 1 (epidemic threshold).
- Highlight the connection: when R_eff crosses 1 downward, it corresponds to the peak of I(t).

#### Scatter Plot — phase portrait
- **Scatter plot**: S vs I (phase plane).
- Show trajectory from the initial point (nearly all S, I ≈ 0) through the peak and back to I = 0.
- Peak point marked with a distinct marker.

### Interactivity

- **Toggle SIR / SEIR** — switch between models.
- **Slider R₀** (0.5–8.0, step 0.1) — the main intuitive parameter; β is recalculated as β = R₀ · γ.
- **Slider γ** (0.01–0.5) — recovery rate (1/γ = average infectious period in days).
- **Slider σ** (0.05–1.0, SEIR only) — rate of symptom onset (1/σ = incubation period in days).
- **Slider "Vaccination"** (0–90%) — initial fraction of immune individuals (R₀ increases, S₀ decreases).
- **Summary panel**: peak infection (day, count), final % recovered, herd immunity threshold (1 − 1/R₀).
- **Disease presets**: buttons "Influenza", "COVID-19", "Measles" that auto-set parameters.

### Tooltips (mandatory, meaningful)

- **R₀**: "Basic reproduction number — average number of people one infected person infects in a fully susceptible population. Influenza ≈ 1.5, COVID-19 ≈ 2.5–3.5, Measles ≈ 12–18."
- **γ (recovery rate)**: "1/γ is the average duration of the infectious period in days. For example, γ = 0.1 means a person is infectious for ~10 days on average."
- **σ (SEIR only)**: "1/σ is the average incubation period in days. During this time the person is infected (E) but not yet infectious to others."
- **Vaccination slider**: "Initial fraction of immune individuals. Herd immunity threshold = 1 − 1/R₀. For R₀ = 3, at least 67% must be vaccinated to prevent an outbreak."
- **Peak of I(t)**: "The infection peak occurs when R_eff = R₀·S/N crosses 1 from above — from that moment each infected person transmits to fewer than one other person on average."
- **Phase portrait S vs I**: "The trajectory always moves leftward (S only decreases in SIR) and ends on the S-axis at I = 0. The final value of S is determined by a transcendental equation."
- **R_eff(t)**: "Effective reproduction number in real time. When it drops below 1, the epidemic is declining."

### Meaningful UI Labels (mandatory)

All headings, axis labels, panel titles, slider labels, and legend entries must be **informative and self-explanatory** — not raw variable names, but phrases that let the user immediately understand what they see and why.

Examples of correct naming:
- Main chart title: **"Epidemic Dynamics — Population Over Time"** — not "SIR Plot", not "Chart"
- Main chart Y-axis: **"Number of individuals"** — not "N", not "count"
- Main chart X-axis: **"Time (days since first case)"** — not "t", not "days"
- Legend entries: **"S — Susceptible (not yet infected)"**, **"I — Infectious (can spread)"**, **"R — Recovered (immune)"**, **"E — Exposed (incubating, SEIR only)"** — not just "S", "I", "R", "E"
- R_eff chart title: **"Effective Reproduction Number R_eff(t)"** — not "R_eff"
- R_eff threshold line label: **"Epidemic threshold (R_eff = 1)"** — not "1.0"
- Phase portrait title: **"Phase Portrait — Susceptible vs Infectious"** — not "S vs I"
- Slider R₀: **"R₀ — basic reproduction number"** — not "R0", not "param"
- Slider γ: **"γ — recovery rate (1/γ = infectious period in days)"** — not "gamma"
- Slider σ: **"σ — incubation rate (1/σ = incubation period in days)"** — not "sigma"
- Vaccination slider: **"Initial vaccination coverage (%)"** — not "vax", not "immune"
- Summary panel heading: **"Epidemic Summary"** — not "Stats", not "Info"
- Peak annotation: **"Peak infection: day 47, 2,340 cases"** — not "max I"
- Preset buttons: **"Influenza (R₀ ≈ 1.5)"**, **"COVID-19 (R₀ ≈ 3.0)"**, **"Measles (R₀ ≈ 15)"** — not "Preset 1"
- Herd immunity line in summary: **"Herd immunity threshold: 67% of population"** — not "HIT = 0.67"

**Principle**: a first-time user opening the simulation should understand every UI element without consulting any documentation.
