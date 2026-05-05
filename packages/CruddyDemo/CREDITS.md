# CruddyDemo — Third-Party Libraries

The `@datagrok/cruddydemo` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

CruddyDemo is a thin demonstration of the `@datagrok-libraries/cruddy` CRUD
app framework. Every runtime dependency is either internal (Datagrok
platform / `@datagrok-libraries/*`) or provided by the platform host as a
webpack external. **No third-party runtime JavaScript is bundled** into the
published artifact, and no copyleft component is bundled.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into CruddyDemo's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published CruddyDemo plugin and impose no obligation on users
of the plugin.
