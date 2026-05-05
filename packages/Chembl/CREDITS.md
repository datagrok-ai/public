# Chembl — Third-Party Libraries

The `@datagrok/chembl` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)).

This package has **no third-party JavaScript dependencies bundled** into its
published artifact. The bulk of Chembl is SQL queries (`queries/*.sql`) and
connection definitions; its TypeScript surface is small and depends only on
the Datagrok JS API and internal `@datagrok-libraries/*` libraries (covered
by the repo-wide MIT umbrella). The cross-package dependency on
`@datagrok/chem` is a `devDependency` only — Chembl never imports from it
at compile time and instead calls `Chem:*` functions through `grok.functions`
at runtime; see [`../Chem/CREDITS.md`](../Chem/CREDITS.md) for Chem's own
third-party attribution.

No copyleft (GPL/LGPL/MPL) component is bundled into the published Chembl
plugin.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

There are no webpack externals declared in Chembl's `dependencies` (the
package does not import `cash-dom`, `dayjs`, `rxjs`, `wu`, or `openchemlib`
directly).

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published Chembl plugin and impose no
obligation on users of the plugin. `file-loader`, `ts-loader`, and
`typescript` are build-time webpack tooling.
