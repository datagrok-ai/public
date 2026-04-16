---
name: create-workflow
description: Create a Datagrok Compute2 workflow (pipeline configuration with steps, links, and actions)
when-to-use: When the user asks to create, design, or scaffold a Compute2 workflow or pipeline configuration. NOT for general scripting, viewers, or UI work.
context: fork
effort: high
argument-hint: "[workflow description]"
---

# Create a Compute2 workflow

Help the user create a Datagrok workflow — a `PipelineConfiguration` object that wires
multiple scripts together with reactive data links, validators, and metadata handlers.

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

1. Read `{SKILL}/references/configuration-reference.md`.
2. Read `{SKILL}/references/link-types-reference.md`.
3. Choose workflow type:
   - `static` — fixed sequence of steps
   - `dynamic` (or `parallel`/`sequential` aliases) — user can add/remove steps at runtime
4. Design the `PipelineConfiguration` object:
   - Define steps with `id` and `nqName`
   - Define data links with `from`/`to` using LQL paths
   - Define validators if cross-field validation is needed
   - Define meta links if UI customization is needed (hiding inputs, setting viewer config)
   - Define actions if user-triggered operations are needed
5. If any link uses tag selectors, template queries, or relative references, read `{SKILL}/references/links-spec-reference.md`.
6. Present the configuration skeleton for approval. Do not implement handlers yet.

### Phase 3: Implement

1. Create the provider function in the package:
   ```typescript
   //name: MyWorkflow
   //description: Description of the workflow
   //output: object config
   //editor: Compute2:TreeWizardEditor
   export function myWorkflow(): PipelineConfiguration {
     return { /* approved configuration */ };
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

Present the review findings. Fix any issues before proceeding.

### Phase 5: Build and verify

1. `grok check --soft` — verify function signatures
2. `webpack` or `npm run build` — build the package
3. `grok publish` — publish to the server
4. Inform the user how to open the workflow in the browser

## Behavior

- **Do not invent scripts.** Only reference functions that exist on the server or in the package.
- **Present config for approval** before writing handler code. Handlers are the expensive part.
- **Use simple LQL paths** unless the user needs dynamic matching. Prefer `in1:step1/a` over complex selectors.
- **Keep handlers pure.** Handlers should transform data, not perform side effects. Use actions for user-triggered operations.
- **One link per output.** Each input of a downstream node should receive data from at most one data link. Use validators or meta links for additional concerns.
- **Use the `/ui` skill** if the workflow needs custom viewer hooks or complex UI layout.
