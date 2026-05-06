# PhyloTreeViewer — Third-Party Libraries

The `@datagrok/phylo-tree-viewer` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is bundled
into the published PhyloTreeViewer plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### deck.gl (8.8.12) and @deck.gl/core (8.8.12)

WebGL2-powered framework for visual exploratory data analysis of large
datasets — provides the rendering substrate for the phylogenetic tree
viewer.

- Upstream: https://github.com/visgl/deck.gl
- License: **MIT**

```
Copyright (c) 2020 Urban Computing Foundation

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

### @phylocanvas/phylocanvas.gl (1.44.0)

WebGL-based phylogenetic tree drawing library — the core engine driving the
tree viewer.

- Upstream: https://gitlab.com/cgps/phylocanvas/-/tree/master/packages/phylocanvas.gl
- License: **MIT**

```
MIT License

Copyright (c) 2021 Centre for Genomic Pathogen Surveillance (CGPS)

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

### mjolnir.js (2.7.x)

Pointer / gesture event manager (drag, pinch, hover) used by the deck.gl
rendering pipeline.

- Upstream: https://github.com/uber-web/mjolnir.js
- License: **MIT**

```
MIT License

Copyright (c) 2017 Uber Technologies, Inc.

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

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into PhyloTreeViewer's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| wu.js     | 2.x     | MIT        | https://github.com/fitzgen/wu.js        |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI, the TypeScript / webpack toolchain, and `puppeteer`
(MIT). The peer/devDependency on `@datagrok/dendrogram` is MIT — covered by
the repo-wide `LICENSE.md`. These are **not** redistributed as part of the
published PhyloTreeViewer plugin and impose no obligation on users of the
plugin.
