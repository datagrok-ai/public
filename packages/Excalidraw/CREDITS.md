# Excalidraw — Third-Party Libraries

The `@datagrok/excalidraw` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is bundled
into the published plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### Excalidraw (@excalidraw/excalidraw, 0.18.0)

Virtual whiteboard / sketching React component used to embed Excalidraw inside
Datagrok. The npm package ships only a `package.json` declaring `"license":
"MIT"`; the canonical MIT license text below is reproduced from the upstream
[`excalidraw/excalidraw`](https://github.com/excalidraw/excalidraw) repository.

The published `dist/` may transitively bundle additional Excalidraw subpackages
(`@excalidraw/laser-pointer`, `@excalidraw/markdown-to-text`,
`@excalidraw/mermaid-to-excalidraw`, `@excalidraw/random-username`); each ships
its own MIT `LICENSE` file in `node_modules/` and carries the same notice.

- Upstream: https://github.com/excalidraw/excalidraw
- Excalidraw application: https://excalidraw.com/
- License: **MIT**

```
MIT License

Copyright (c) 2020 Excalidraw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

**Embedded fonts.** Excalidraw bundles its own font files (Cascadia, Virgil,
Lilita One, Comic Shanns, Nunito, etc.) under `dist/prod/fonts/`. Each font
carries its own SIL Open Font License (OFL) or other permissive license; the
license texts are shipped alongside the font binaries inside the Excalidraw
distribution. By default the Datagrok bundle relies on the Excalidraw CDN
(`window.EXCALIDRAW_ASSET_PATH` unset) and does not redistribute the font
binaries directly.

### roughjs (4.6.x)

Hand-drawn-style rendering library used by Excalidraw to produce its
distinctive sketchy shapes.

- Upstream: https://github.com/rough-stuff/rough — https://roughjs.com/
- License: **MIT**

```
MIT License

Copyright (c) 2019 Preet Shihn

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### React / React DOM (transitive — peer of @excalidraw/excalidraw)

Excalidraw is a React component, so React and React-DOM are pulled in as peer
dependencies. License: **MIT** — *Copyright (c) Meta Platforms, Inc. and
affiliates.* See `node_modules/react/LICENSE` for the full text.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Excalidraw's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Excalidraw plugin and impose no obligation on users of
the plugin.
