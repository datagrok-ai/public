# Compute — Third-Party Libraries

The `@datagrok/compute` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

Compute is the legacy modeling UI. Every runtime dependency is either internal
(Datagrok platform / `@datagrok-libraries/*`, including `compute-utils` which
holds the bulk of the modeling code) or provided by the platform host as a
webpack external. **No third-party runtime JavaScript is bundled** into the
published artifact, and no copyleft component is bundled.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Compute's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| ExcelJS                          | 4.x     | MIT          | https://github.com/exceljs/exceljs                |

`exceljs` and `html2canvas` are also externalised in `webpack.config.js` as a
defensive measure even though they are not currently imported by Compute's TS
sources; they are kept on the externals list so any feature that imports them
will pick up the platform copy.

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Compute plugin and impose no obligation on users of
the plugin.
