---
name: create-interactive-scientific-application
description: >
  Build interactive scientific applications on the Datagrok platform using a
  structured two-phase workflow: first create a detailed specification from a
  template, then implement following the architecture guide. Use this skill
  whenever asked to build simulations, interactive models, scientific tools,
  ODE/PDE solvers, data analysis apps, or any Datagrok package application.
  Also use when the user mentions "Datagrok app", "interactive app",
  "scientific application", "spec template", "specification template",
  or asks to create a Datagrok package function with UI and computations.
---

# Datagrok Interactive Scientific Application Development

This skill implements a two-phase workflow for building production-quality
interactive scientific applications on the Datagrok platform.

**Phase 1** — Build a detailed specification from a template.
**Phase 2** — Implement the application following the architecture guide.

Never skip Phase 1. Never start coding before the specification is approved.

**Base path for all references in this skill:**
`.claude/skills/create-interactive-scientific-application/` (relative to repo root).
All file paths below use `{SKILL}/` as shorthand for this base path.

---

## Phase 1: Specification

### Step 1: Understand the request

Analyze the user's request and identify:
- The scientific domain and problem being solved.
- Core computational model (formulas, equations, algorithms).
- What the user wants to see and interact with.
- Complexity level: simple (single reactive task, no workers) vs. complex
  (multiple tasks, workers, secondary pipelines).

Determine the target package:
- If the user specifies a package name — use it.
- If the current working directory is inside an existing package
  (has `package.json` with `datagrok-api` dependency) — add the
  application to that package.
- Otherwise — ask the user whether to create a new package or add
  to an existing one.

If the request is vague (e.g., "build a Lotka-Volterra app"), ask clarifying
questions before proceeding — but keep it to 2–3 focused questions, not an
interrogation.

### Step 2: Read reference materials

1. Read the spec template in full: `{SKILL}/references/spec-template.md`
2. Check what completed examples exist:
   `{SKILL}/references/examples/` — list available subdirectories.

Do NOT read the guide or examples in full upfront. Read them in stages,
aligned with the sections you are currently filling (Step 3) and the
approval flow (Step 4):

- When filling sections 1.0–1.5 → read guide sections 1–1.5, then
  example sections 1–1.5.
- When filling sections 2–5 → read guide sections 2–5, then
  example sections 2–5.
- When filling sections 6–15 → read guide sections 6–15, then
  example sections 6–15.

File paths:
- Guide: `{SKILL}/references/guide.md`
- Example spec (has a TOC with line ranges for partial reading):
  look inside `{SKILL}/references/examples/` subdirectories.

### Step 3: Fill in the specification

Fill in every section of the template. Follow these rules:

#### Critical sections — fill with maximum detail:

- **1.0 General Information** — application name, package, entry function, description.
- **1.1 Core** — task list, computation formulas (Level 1 is MANDATORY before
  approval), input/output parameters, implementation approach.
- **3. Controls** — every control must have: ID, label, type, data type,
  default value, min, max, format, tooltip. No placeholders.
- **4. Display Elements** — every viewer/element with ID, type, associated data.
- **7. Validation** — concrete rules with conditions and error messages.
- **15. Testing** — test categories, reference examples, expected coverage.

#### Sections that may be N/A for simple apps:

For applications with a single reactive task, no workers, and no secondary
pipelines, the following sections can be marked N/A with a one-line explanation:
- 1.2 Ports (application-level: Progress, Cancellation)
- 8.2 Secondary Pipelines
- 8.4 Computation Blocking
- 12.2 Worker Termination

Never mark a section N/A without explanation. If unsure whether a section
applies — include it.

#### Computation formulas (Section 1.1) — special attention:

This is the most important section. Level 1 must contain:
- All variables with meaning, units, and valid domains.
- All equations/relationships connecting inputs to outputs — unambiguously,
  so that another developer could implement from this description alone.
- Output properties (invariants): bounds, monotonicity, conservation laws,
  limiting cases.
- At least one reference example per computational path: concrete inputs →
  expected output with source (manual calculation / literature).

Level 2 (full formalization) can be marked as "to be developed incrementally"
for the first iteration, but note what it will eventually contain.

#### Control IDs — naming convention:

Use a consistent prefix: `ctrl_` for inputs (e.g., `ctrl_alpha`, `ctrl_x0`),
`btn_` for buttons (e.g., `btn_optimize`, `btn_reset`), `view_` for viewers
(e.g., `view_line_chart`, `view_phase_plot`).

### Step 4: Present for approval

Present the specification to the user in stages:
1. **First:** Sections 1.0–1.5 (architecture, computation model, ports, adapters,
   coordinator, independence principle).
   Ask for feedback before continuing — this is the foundation.
2. **Then:** Sections 2–5 (main view, controls, display, layout).
3. **Then:** Sections 6–15 (feedback, validation, pipeline, reactivity,
   data lifecycle, errors, resources, closure, UX, testing).

At each stage, ask: "Does this look correct? Anything to change?"

When the user approves the full specification — confirm explicitly:
"Specification approved. Moving to implementation."

### Iteration rules

- If the user requests changes — update the spec and re-present the changed
  sections.
- If you discover ambiguities or inconsistencies while filling in later
  sections — go back and fix earlier sections, noting the changes.
- The specification is the single source of truth for Phase 2.

---

## Phase 2: Implementation

### Step 1: Read the architecture guide

Read `{SKILL}/references/guide.md` in full if not already loaded.

If the guide references additional reference documents that are relevant
to this application, check `{SKILL}/references/reference/` for them.
The guide may reference files like:
- `COMPUTATION-PATTERNS.md` — patterns for raw data and null handling.
- `ARRAY-OPERATIONS.md` — efficient typed array operations.
- `WORKER-GUIDE.md` — worker-utils infrastructure and lifecycle.
- `PARALLEL-EXECUTION.md` — distribution across multiple workers.

Additionally, always read these before implementing UI and writing code:
- `datagrok-api-reference.md` — catalog of Datagrok inputs, viewers,
  buttons, tooltips, dialogs, layouts, notifications, subscriptions.
- `datagrok-coding-conventions.md` — file structure, naming, formatting,
  error handling, testing conventions.

Read those that are relevant to the current application.

### Step 2: Implement in order

Follow this implementation sequence strictly:

1. **Project structure** — create files and directories per the guide.
2. **Model types** — interfaces for input/output ports (Section 1.2 of spec).
3. **Core computation** — implement each task in isolation, no UI imports
   (Section 1.1 of spec). The core must NOT import `datagrok-api` or `ui`.
4. **Validation** — implement `validate()` returning `Map<InputId, string>`
   (Section 7 of spec).
5. **Workers** (if applicable) — worker files in `src/<app-name>/workers/`
   (Section 1.3 of spec). Each worker must be added as a separate entry
   point in `webpack.config.js`:
   ```
   '<worker-name>': {filename: '<worker-name>.js', import: './src/<app-name>/workers/<worker-name>.ts'}
   ```
   Workers are loaded at runtime via `new Worker(_package.webRoot + 'dist/<worker-name>.js')`.
6. **UI controls** — create inputs with all options from the spec
   (Section 3 of spec).
7. **Display elements** — viewers, custom panels (Section 4 of spec).
8. **CSS** — all styles in `css/<app-name>.css` with app-specific prefix.
   No inline styles. (Section 5.3 of spec).
9. **Coordinator** — connect everything: reactivity, validation triggering,
   computation, result display, resource lifecycle (Section 1.4 of spec).
10. **Tests** — implement test categories from Section 15 of spec.

### Step 3: Implementation rules

Follow the architecture guide strictly. Key principles:

- **Hexagonal architecture:** Core ↔ Ports ↔ Adapters ↔ Coordinator.
- **Independence principle:** UI behavior does not depend on computation.
  The core receives a ready, validated parameter set.
- **CSS isolation:** all classes use `<app-name>-` prefix. No generic names.
- **Computation blocking:** use `computationsBlocked` flag for batch updates.
- **Resource cleanup:** collect all subscriptions in `subs[]`, terminate
  workers on close, cancel pending operations.
- **Tooltips:** every control and action button must have a tooltip.
- **Validators:** add Datagrok validators to every input via `addValidators()`.
- **No artificial delays:** the example `equilibrium-worker.ts` contains a
  `busy wait` loop for demo purposes — never reproduce this pattern.
- **Domain names in core:** the example `core.ts` uses UI-prefixed
  identifiers (`ctrl_alpha`, `ctrl_beta`) in `InputId` for brevity.
  In production code, the core should use domain names (`alpha`, `beta`)
  and the coordinator should map them to control IDs. Similarly, `RANGES`
  (slider min/max) belongs to the UI adapter, not the core — the core
  should only define mathematical domains (e.g., `> 0`).

### When the spec is insufficient

If during implementation you discover that the specification is incomplete,
ambiguous, or inconsistent:
- **STOP coding.**
- Tell the user exactly what is missing or unclear.
- Propose a fix.
- Wait for approval before continuing.

Do NOT improvise or fill gaps silently.

---

## Quality checklist (before presenting the final result)

- [ ] Every computation matches the formulas in the spec exactly.
- [ ] All controls have correct defaults, ranges, formats, and tooltips.
- [ ] Validation covers all rules from the spec with correct error messages.
- [ ] CSS file exists with app-prefixed classes; no inline styles.
- [ ] All event subscriptions collected in `subs[]` array.
- [ ] All workers terminated on close.
- [ ] `onViewRemoved` handler cleans up all resources.
- [ ] Tests cover: validation (boundary + invalid), formula verification,
      output properties, numerical method (if applicable).
