# PowerGrid — Third-Party Libraries

The `@datagrok/power-grid` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

This package has no third-party runtime JavaScript dependencies bundled into
its published artifact — all runtime libraries are either Datagrok-internal
(`@datagrok-libraries/*`, `datagrok-api`) and covered by the repo-wide MIT
`LICENSE.md`, or webpack externals provided by the platform host (Section 2).

---

## 1. Bundled in the published artifact (`dist/`)

None.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into PowerGrid's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |
| wu.js     | 2.x     | MIT     | https://github.com/fitzgen/wu.js        |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack / Babel toolchain. These
are **not** redistributed as part of the published PowerGrid plugin and impose
no obligation on users of the plugin.
