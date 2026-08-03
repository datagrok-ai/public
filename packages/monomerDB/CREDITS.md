# monomerDB — Third-Party Libraries

The `@datagrok/monomerdb` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

monomerDB has no third-party runtime dependencies of its own — every external
library it uses is either an internal Datagrok dependency (`@datagrok-libraries/*`,
covered by the repo-wide MIT) or a webpack external provided by the platform.

---

## 1. Bundled in the published artifact (`dist/`)

None.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into monomerDB's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Docker container images (`dockerfiles/`)

None. monomerDB does not ship a Docker image.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published monomerDB plugin and impose no obligation on users
of the plugin.
