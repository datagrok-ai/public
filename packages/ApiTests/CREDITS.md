# ApiTests — Third-Party Libraries

The `@datagrok/api-tests` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)).

This package has **no third-party JavaScript dependencies bundled** into its
published artifact. All non-Datagrok runtime dependencies declared in
`package.json` (`cash-dom`, `dayjs`, `rxjs`, `wu`) are webpack externals —
they are not bundled into ApiTests' `dist/`, but are provided once by the
Datagrok platform host and shared across all packages.

No copyleft (GPL/LGPL/MPL) component is bundled into the published ApiTests
plugin.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |
| wu.js     | 2.x     | MIT        | https://github.com/fitzgen/wu.js        |

---

## 2. Docker container (`dockerfiles/apitests-docker-test1`) — Test fixture

ApiTests ships a tiny test-fixture Docker image (`dockerfiles/apitests-docker-test1/Dockerfile`)
used by the test suite to exercise the `grok.dapi.docker` API. It is built
from `python:3.8-alpine3.18` and installs:

| Component | License     |
|-----------|-------------|
| Python    | PSF License |
| Flask     | BSD-3-Clause |

This image is a test fixture only and is not published to end users.

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published ApiTests plugin and impose no
obligation on users of the plugin.

The `yargs` package is declared as an `optionalDependencies` entry; it is
used by the Node-based offline test runner (`start-node`/`stress-node`), not
by the webpack bundle.
