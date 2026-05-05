# GIS — Third-Party Libraries

The `@datagrok/gis` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (BSD-2-Clause, MIT). No copyleft (GPL/LGPL/MPL) component
is bundled into the published plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### OpenLayers (ol, 6.15.x)

The map rendering engine that drives all GIS map viewers in the package.

- Upstream: https://openlayers.org/ — https://github.com/openlayers/openlayers
- License: **BSD-2-Clause**

```
BSD 2-Clause License

Copyright 2005-present, OpenLayers Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### JSZip (3.10.x)

ZIP file reader/writer used to import zipped Shapefiles and KMZ archives.

- Upstream: https://stuk.github.io/jszip/
- License: **MIT** (chosen from upstream's `MIT OR GPL-3.0-or-later` dual license)

```
Copyright (c) 2009-2016 Stuart Knightley, David Duponchel, Franz Buchinger, António Afonso

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

### xhr2 (0.2.x)

XMLHttpRequest implementation used in build / test scripts that load remote
tiles outside a browser.

- Upstream: https://github.com/pwnall/node-xhr2
- License: **MIT**

```
Copyright (c) 2013 Victor Costan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into GIS's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

OpenLayers map viewers fetch raster/vector tiles from third-party tile
providers (OpenStreetMap, Stadia Maps, Carto, Esri, etc.) at runtime,
according to the basemap configuration chosen by the user. These tiles are
**not** bundled into GIS's `dist/`. Each provider's terms of use, attribution
and usage limits apply to the deployed Datagrok instance.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published GIS plugin and impose no obligation on users of the
plugin.
