---
paths:
  - "**/webpack.config.*"
---

## Webpack Configuration

Packages typically have two entry points: `package.ts` (main) and `package-test.ts` (tests).

These must be marked as externals (provided by platform at runtime):
- `datagrok-api` (including subpaths `/grok`, `/ui`, `/dg`)
- `rxjs`, `rxjs/operators`
- `cash-dom`, `dayjs`, `wu`
- `openchemlib/full.js`

Use `grok check` to validate that imports match webpack externals.
