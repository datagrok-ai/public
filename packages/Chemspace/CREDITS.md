# Chemspace — Third-Party Libraries

The `@datagrok/chemspace` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This package integrates with the Chemspace REST API (see [`CLAUDE.md`](CLAUDE.md)).
It bundles **no third-party runtime JavaScript** into the published artifact:
every runtime dependency is either internal (the Datagrok platform /
`@datagrok-libraries/*`) or provided by the platform host as a webpack external.
No copyleft component is bundled.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Chemspace's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |

Only the three `datagrok-api` subpaths (`/dg`, `/grok`, `/ui`) are externalised
in this package's `webpack.config.js`; the other repo-wide externals
(`rxjs`, `cash-dom`, `dayjs`, `wu`, `openchemlib/full.js`) are not imported by
this package.

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Chemspace plugin and impose no obligation on users of
the plugin.
