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

## Development

```bash
npm install
npm run build      # tsc, compiles in place under src/
npm run gallery    # serves the standalone gallery at /gallery/
```
