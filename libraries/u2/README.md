# @datagrok-libraries/u2

Next-gen Datagrok UI library: design tokens, a headless signal-based component core, and a
standalone component gallery. Compact, IDE-class visual language by default.

Architecture (see the design docs in the core repo, `core/docs/features/ui2/`):

- **L0 tokens** — CSS custom properties, single source of truth for both old CSS and u2.
- **L1 headless core** — TypeScript + signals; state, behavior, a11y, validation.
- **L2 spec & registry** — JSON-serializable UI, machine-readable component manifest.
- **L3 delivery** — fluent `u2.*` API, custom-element tags, generated framework wrappers.

L0–L2 are **platform-free**: no `datagrok-api` imports (enforced by eslint). All platform
integration — `DG.Widget` lifecycle, `InputBase` conformance, Property/semType-driven
inputs — lives in the `@datagrok-libraries/u2/dg` entry point.

`vendor/` holds frozen copies of permissively-licensed upstream code; see
[vendor/LICENSES.md](vendor/LICENSES.md).

## Layout

```
src/core/          signals, Scope, Component/Control, Input base, elements, async primitives
src/components/    platform-free controls, by role — the same vocabulary as the registry categories
  inputs/          everything with a value the form and the Dart bridge read and write (+ combobox, typeahead)
  forms/           form, property-grid — edit an object through rows of inputs
  containers/      splitter, tabs, accordion, dialog, wizard — arrange children, carry no value
  actions/         buttons, button-group, row/context actions
  navigation/      menu, menu-bar, toolbar, breadcrumbs
  collections/     virtualized list, tree, table
  display/         icon, badge, progress, async-view
src/spec/          dg-ui/1 document: renderer, registry, patch engine, binding paths
src/sources/       the data sources a spec can declare (platform-free, behind `backends`)
src/dg/            the platform layer — the only code that imports datagrok-api; `dg/designer/` is the visual designer
css/               one sheet per component, flat — imported by path from packages
```

The root `src/index.ts` is the single public surface; import from it (or from `src/dg/index.ts`),
never from a component's file path.

## Development

```bash
npm install
npm run build      # tsc, compiles in place under src/
npm run gallery    # serves the standalone gallery at /gallery/
```
