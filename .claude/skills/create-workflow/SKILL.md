---
name: create-workflow
description: Create a Datagrok Compute2 workflow (pipeline configuration with steps, links, and actions)
when-to-use: When the user asks to create, design, or scaffold a Compute2 workflow or pipeline configuration. NOT for general scripting, viewers, or UI work.
context: fork
effort: high
argument-hint: "[workflow description]"
disable-model-invocation: true
---

# Create a Compute2 workflow

The user wants a `PipelineConfiguration` — a Datagrok Compute2 workflow that wires
multiple scripts together with reactive data links, validators, and metadata handlers.

The authoritative reference lives in the docs. This skill is a working procedure;
when you need to know _what something is_ or _how a field behaves_, read the docs.

## Reading order

Read these before you write any configuration:

- `help/compute/workflows/overview.mdx` — terms (node, link, controller, action, FuncCall, nqName, RichFunctionView).
- `help/compute/workflows/configuration.mdx` — every field of `PipelineConfiguration`,
  the workflow types, states, custom exports, and the
  [constraints / review checklist](#) (consult the section before publishing).
- `help/compute/workflows/link-types.mdx` — link/action types, controller methods, handler signatures.

Read on demand:

- `help/compute/workflows/links-spec.mdx` — Link Query Language grammar. Needed only when the
  workflow uses tag selectors, template queries, or relative (`base` / `@base`) refs.
- `help/compute/workflows/examples.mdx` — Wine Quality walkthrough end-to-end.
- `help/compute/workflows/code-usage.mdx` — only if the workflow will be launched
  programmatically.

Reference examples (also linked from `examples.mdx`):

| File | Use when |
|------|----------|
| `examples/minimal-static.ts` | Fixed sequence of scripts, no links, user fills inputs manually. |
| `examples/dynamic-with-links.ts` | User can add/remove steps; outputs propagate to all downstream instances. |
| `examples/validators-and-meta.ts` | Cross-field validation, conditional input visibility, user-triggered actions. |

Setup (install + dayjs/timezone imports) is covered in `examples.mdx#dependencies`.

## Instructions

### Phase 1: Understand the requirements

1. Ask the user:
   - What scripts (Datagrok functions) should the workflow connect?
   - What data flows between them? (which outputs feed which inputs)
   - Is the set of steps fixed (static) or user-configurable (dynamic)?
   - Are there validation rules? (required fields, value ranges, cross-field checks)
   - Should any inputs be visually customized? (hidden, readonly, dropdowns)
2. Verify that every referenced script exists. A script may be deployed, scaffolded
   locally but not yet published, or only an idea in the user's head. Check each
   location and stop searching once you find a match:
   - **Deployed on the server**: `grok s functions list --filter "<nqName>"`.
     A non-empty result means the script is live and the `nqName` is correct.
   - **Local `package.ts`**: grep for `//name:\s*<FunctionName>` in `src/package.ts`
     (and any `src/package-*.ts` entries). Each annotated export becomes a function
     with `nqName: <PackageName>:<FunctionName>` once published.
   - **Local `scripts/` directory**: grep for `^#name:\s*<scriptName>` in
     `scripts/**/*.{py,r,js,jl,m,sql}`. Each `#name`-annotated file becomes a
     function with `nqName: <PackageName>:<scriptName>` once published.

   If a script is found locally but not on the server, note it as "scaffolded — will be
   published with this workflow". If a script is missing in all three places, ask the
   user whether to scaffold it (and follow the appropriate skill: see `/init` or the
   scripting docs) or to drop it from the workflow.
3. Present a plain-language summary of the workflow for approval before coding. Mark
   each step as **deployed**, **scaffolded**, or **to be created** so the user can see
   the integration surface at a glance.

### Phase 2: Design the configuration

1. Choose the workflow type. See `configuration.mdx` for the discriminated union of
   `static` / `dynamic` / `action` / `ref`.
2. Sketch the `PipelineConfiguration` object: steps with `id` and `nqName`; data links
   with `from`/`to`; validators and meta links if needed; actions for user-triggered
   operations.
3. If any link uses tag selectors, template queries, or relative references, consult
   `links-spec.mdx`.
4. Present the configuration skeleton for approval. Do not implement handlers yet.

### Phase 3: Implement

1. Create the provider function in the package:
   ```typescript
   import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';

   //name: MyWorkflow
   //description: Description of the workflow
   //tags: model
   //editor: Compute2:TreeWizardEditor
   //input: object params
   //output: object result
   export function myWorkflow(): PipelineConfiguration {
     return {
       id: 'my-workflow',
       nqName: 'MyPackage:MyWorkflow',
       version: '1.0',
       /* approved configuration */
     };
   }
   ```
2. Implement link handlers using the controller methods documented in `link-types.mdx`.
3. Register the function in `package.ts` if not already there.
4. Run `grok api` to regenerate wrappers.

### Phase 4: Review

Validate the configuration against the
[constraints and review checklist](../../../help/compute/workflows/configuration.mdx#constraints-and-review-checklist).
Spawn a sub-agent for an independent pass if the configuration is non-trivial. Fix
anything that fails before proceeding.

### Phase 5: Build and verify

1. `grok check --soft` — verify function signatures.
2. `webpack` or `npm run build` — build the package.
3. `grok publish --release` — `--release` is mandatory; debug-mode packages are only
   visible to the publishing user.
4. Tell the user where to open the workflow: **Apps → Compute → ModelHub**.

## Behavior

- **Do not invent scripts.** Only reference functions that exist on the server or in the package.
- **Present config for approval** before writing handler code. Handlers are the expensive part.
- **Use simple LQL paths** unless the user needs dynamic matching. Prefer `in1:step1/a`
  over complex selectors.
- **Keep handlers pure.** Handlers should transform data, not perform side effects.
  Use actions for user-triggered operations.
- **One data link per script input.** Each input of a downstream node should receive
  data from at most one data link. Use validators or meta links for additional concerns.
- **Import from `@datagrok-libraries/compute-api`.** This is the public API. Do not
  import from `@datagrok-libraries/compute-utils` directly — those are internal paths.
- **Always publish with `--release`.**
- **Use the `/ui` skill** only in two cases: (1) an action needs to show custom inputs
  inline (e.g. a confirmation form with extra fields), or (2) a script's output viewer
  needs tweaks applied through its `DG.Viewer` JS API inside a
  [`viewersHook`](../../../help/compute/workflows/configuration.mdx#viewershook-script-node).
  Workflow scaffolding, links, validators, and meta-driven UI changes do not need `/ui`.
