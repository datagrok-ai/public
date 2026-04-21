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

Help the user create a Datagrok workflow — a `PipelineConfiguration` object that wires
multiple scripts together with reactive data links, validators, and metadata handlers.

## Prerequisites

The package must have `@datagrok-libraries/compute-api` installed. This is the public API
for workflow development — it re-exports all needed types (`PipelineConfiguration`, `Handler`,
`Validator`, `MetaHandler`, controller interfaces, etc.).

```bash
npm i @datagrok-libraries/compute-api
```

The provider file must also import `dayjs` with UTC/timezone plugins — required by
TreeWizardEditor at runtime:

```typescript
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import timezone from 'dayjs/plugin/timezone';
dayjs.extend(utc);
dayjs.extend(timezone);
```

## References

**Always read before implementing:**

| File | Purpose |
|------|---------|
| `help/compute/workflows/configuration.mdx` | All configuration options, workflow types, states, exports, edge cases |
| `help/compute/workflows/link-types.mdx` | Link/action types, controller methods, handler signatures |

**Read when the workflow uses complex link queries (tags, templates, relative refs):**

| File | Purpose |
|------|---------|
| `help/compute/workflows/links-spec.mdx` | Link Query Language (LQL) specification |

**Read for additional context as needed:**

| File | Purpose |
|------|---------|
| `help/compute/workflows/overview.mdx` | Core terms and concepts (node, link, controller, action, FuncCall) |
| `help/compute/workflows/examples.mdx` | Full tutorial: Wine Quality dataset workflow with 5 steps, data links, visualization |
| `help/compute/workflows/code-usage.mdx` | Starting workflows programmatically from code |

**Read the example that matches the target complexity:**

| File | Complexity |
|------|-----------|
| `.claude/skills/create-workflow/examples/minimal-static.ts` | Static pipeline, no links |
| `.claude/skills/create-workflow/examples/dynamic-with-links.ts` | Dynamic pipeline with data links and restrictions |
| `.claude/skills/create-workflow/examples/validators-and-meta.ts` | Validators, meta handlers, and actions |

## Instructions

### Phase 1: Understand the requirements

1. Ask the user:
   - What scripts (Datagrok functions) should the workflow connect?
   - What data flows between them? (which outputs feed which inputs)
   - Is the set of steps fixed (static) or user-configurable (dynamic)?
   - Are there validation rules? (required fields, value ranges, cross-field checks)
   - Should any inputs be visually customized? (hidden, readonly, dropdowns)
2. Verify that the referenced scripts exist: `grok s functions list --filter "<nqName>"`
3. Present a plain-language summary of the workflow for approval before coding.

### Phase 2: Design the configuration

1. Read `help/compute/workflows/configuration.mdx`.
2. Read `help/compute/workflows/link-types.mdx`.
3. Choose workflow type:
   - `static` — fixed sequence of steps
   - `dynamic` — user can add/remove steps at runtime
   - `action` — lightweight placeholder for displaying actions via visibleOn
4. Design the `PipelineConfiguration` object:
   - Define steps with `id` and `nqName`
   - Define data links with `from`/`to` using LQL paths
   - Define validators if cross-field validation is needed
   - Define meta links if UI customization is needed (hiding inputs, setting viewer config)
   - Define actions if user-triggered operations are needed
5. If any link uses tag selectors, template queries, or relative references, read `help/compute/workflows/links-spec.mdx`.
6. Present the configuration skeleton for approval. Do not implement handlers yet.

### Phase 3: Implement

1. Create the provider function in the package:
   ```typescript
   import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';

   //name: MyWorkflow
   //description: Description of the workflow
   //tags: model
   //output: object config
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
2. Implement link handlers using the controller methods from the reference.
3. Register the function in `package.ts` if not already there.
4. Run `grok api` to regenerate wrappers.

### Phase 4: Review

Spawn a sub-agent to review the configuration. The reviewer must check:

- [ ] **No cycles:** data links only reference the same or forward nodes in DFS order
- [ ] **LQL paths resolve:** every `from`/`to` path segment matches a step `id` in the config
- [ ] **Controller methods match link type:** `setAll` only in data/action handlers, `setValidation` only in validators, `setViewMeta` only in meta handlers
- [ ] **Restrictions consistent:** `defaultRestrictions` types match the intended behavior (`restricted` = editable with warning, `disabled` = locked, `info` = informational)
- [ ] **nodePriority set** when multiple links write to the same node and order matters
- [ ] **No duplicate I/O deps:** at most one data link writes to each input of a node
- [ ] **Import path:** types imported from `@datagrok-libraries/compute-api`, not from `compute-utils`

Present the review findings. Fix any issues before proceeding.

### Phase 5: Build and verify

1. `grok check --soft` — verify function signatures
2. `webpack` or `npm run build` — build the package
3. `grok publish --release` — publish to the server (must use `--release` for workflows)
4. Inform the user how to open the workflow: **Apps -> Compute -> ModelHub**

## Behavior

- **Do not invent scripts.** Only reference functions that exist on the server or in the package.
- **Present config for approval** before writing handler code. Handlers are the expensive part.
- **Use simple LQL paths** unless the user needs dynamic matching. Prefer `in1:step1/a` over complex selectors.
- **Keep handlers pure.** Handlers should transform data, not perform side effects. Use actions for user-triggered operations.
- **One link per output.** Each input of a downstream node should receive data from at most one data link. Use validators or meta links for additional concerns.
- **Import from `@datagrok-libraries/compute-api`.** This is the public API. Do not import from `@datagrok-libraries/compute-utils` directly — those are internal paths.
- **Always publish with `--release`.** Debug-mode packages are only visible to the developer.
- **Use the `/ui` skill** if the workflow needs custom viewer hooks or complex UI layout.
