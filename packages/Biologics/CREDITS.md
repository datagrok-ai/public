# Biologics — Third-Party Libraries

The `@datagrok/biologics` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)).

This package has **no third-party JavaScript dependencies bundled** into its
published artifact. All non-Datagrok runtime dependencies declared in
`package.json` (`cash-dom`, `dayjs`) are webpack externals — they are not
bundled into Biologics' `dist/`, but are provided once by the Datagrok
platform host and shared across all packages.

No copyleft (GPL/LGPL/MPL) component is bundled into the published Biologics
plugin.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published Biologics plugin and impose no
obligation on users of the plugin.
