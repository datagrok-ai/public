# Tableau — Third-Party Libraries

The `@datagrok/tableau` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This package adds import / preview support for Tableau `.twb` workbook files.
It contains no bundled third-party JavaScript runtime libraries beyond what
the Datagrok platform already provides as webpack externals; .twb files are
parsed as plain XML using the browser's built-in DOM parser.

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either listed in Section 2 (provided by the
platform) or are Datagrok-internal libraries covered by the repo-wide
`LICENSE.md`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Tableau's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                          |
|-----------|---------|------------|-----------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs |

---

## 3. Fetched at runtime from third-party services

None. No Tableau SDK or proprietary Tableau JavaScript API is bundled or
loaded at runtime; this plugin only reads local `.twb` XML files.

---

## 4. Docker container images

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Tableau plugin and impose no obligation on users of
the plugin.
