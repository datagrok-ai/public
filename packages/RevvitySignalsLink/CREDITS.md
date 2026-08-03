# RevvitySignalsLink — Third-Party Libraries

The `revvitysignalslink` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This package is a thin REST-API client for the Revvity Signals platform; it
contains no bundled third-party JavaScript runtime libraries beyond what the
Datagrok platform already provides as webpack externals.

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either listed in Section 2 (provided by the
platform) or are Datagrok-internal libraries covered by the repo-wide
`LICENSE.md`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into RevvitySignalsLink's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party services (not bundled)

The package talks to a customer-configured **Revvity Signals** REST endpoint
through `grok.dapi.fetchProxy`. No Revvity SDK or JS code is bundled or
redistributed by this plugin; access to the Revvity Signals service is
governed entirely by the customer's own Revvity license.

---

## 4. Docker container images

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published RevvitySignalsLink plugin and impose no obligation on
users of the plugin.
