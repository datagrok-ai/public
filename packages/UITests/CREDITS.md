# UITests — Third-Party Libraries

The `@datagrok/ui-tests` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This is an internal automated-UI-tests service package; it contains no bundled
third-party JavaScript runtime libraries beyond what the Datagrok platform
already provides as webpack externals.

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either listed in Section 2 (provided by the
platform) or are Datagrok-internal libraries covered by the repo-wide
`LICENSE.md`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into UITests' `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |
| wu.js     | 2.x     | MIT     | https://github.com/fitzgen/wu.js        |

---

## 3. Fetched at runtime from third-party services

None.

---

## 4. Docker container images

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published UITests plugin and impose no obligation on users of
the plugin.
