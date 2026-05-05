# Notebooks — Third-Party Libraries

The `@datagrok/notebooks` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

The plugin embeds large parts of **JupyterLab** and the **Lumino** widget
toolkit to render Jupyter notebooks inside the Datagrok UI. All bundled
JavaScript dependencies are under permissive licenses (MIT, BSD-3-Clause).
No copyleft (GPL/LGPL/MPL) JavaScript component is bundled into the published
Notebooks plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### JupyterLab packages (2.3.x) — `@jupyterlab/{apputils, codemirror, completer, coreutils, docmanager, docregistry, documentsearch, mathjax2, notebook, rendermime, services, theme-light-extension}`

The frontend pieces of JupyterLab needed to host a notebook view in the
browser. All `@jupyterlab/*` packages share the same upstream license.

- Upstream: https://github.com/jupyterlab/jupyterlab
- Copyright: Project Jupyter Contributors
- License: **BSD-3-Clause**

```
Copyright (c) 2015 Project Jupyter Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

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

The `@jupyterlab/coreutils` tarball includes `semver.py`, which is licensed
under MIT (see the file header).

### Lumino (`@lumino/widgets`, `@lumino/commands`, …)

Phosphor / Lumino is JupyterLab's widget toolkit, also pulled in transitively
by every `@jupyterlab/*` package.

- Upstream: https://github.com/jupyterlab/lumino
- Copyright: Project Jupyter Contributors
- License: **BSD-3-Clause**

> Lumino is distributed under the same BSD-3-Clause text shown above for
> JupyterLab. The full text is available at
> https://github.com/jupyterlab/lumino/blob/main/LICENSE.

### bufferutil (4.x)

WebSocket buffer-manipulation helper (used by JupyterLab kernel transport).

- Upstream: https://github.com/websockets/bufferutil
- License: **MIT**

```
Copyright (c) 2011 Einar Otto Stangvik <einaros@gmail.com>
Copyright (c) 2013 Arnout Kazemier and contributors
Copyright (c) 2016 Luigi Pinca and contributors

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

### utf-8-validate (5.x)

UTF-8 validator companion to `bufferutil`.

- Upstream: https://github.com/websockets/utf-8-validate
- License: **MIT** (analogous notice to `bufferutil` above)

### es6-promise (4.x)

Promise polyfill for older runtimes.

- Upstream: https://github.com/stefanpenner/es6-promise
- License: **MIT**

```
Copyright (c) 2014 Yehuda Katz, Tom Dale, Stefan Penner and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

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

### Other transitive bundled libraries

JupyterLab transitively pulls in dozens of small permissive dependencies
(CodeMirror, MathJax 2, Marked, Sanitize, ajv, etc.). The set has been
audited to be MIT / BSD / Apache-2.0. The `@jupyterlab/mathjax2` package
includes the MathJax 2 distribution, licensed under **Apache-2.0**.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

None — Notebooks does not import any of the standard externals.

---

## 3. Docker container images (`dockerfiles/`)

### `dockerfiles/jupyter-notebook` — Jupyter notebook server with nginx proxy

The image is built `FROM datagrok/jupyter_notebook:bleeding-edge` (a Datagrok-
maintained image) and adds an `nginx` reverse proxy + custom Jupyter config
files.

| Component                   | Source                                        | License                  |
|-----------------------------|-----------------------------------------------|--------------------------|
| Jupyter Notebook / IPython  | https://jupyter.org/                          | BSD-3-Clause             |
| Python                      | https://www.python.org/                       | PSF                      |
| nginx                       | https://nginx.org/                            | BSD-2-Clause             |
| `gettext`                   | https://www.gnu.org/software/gettext/         | GPL-3.0+ (utilities only)|
| Debian base packages        | https://www.debian.org/                       | various                  |

The Debian base image and the `datagrok/jupyter_notebook` upstream image
include GNU/Linux system libraries under copyleft licenses (glibc — LGPL;
coreutils, bash, apt, gettext utilities — GPL). These are standard OS
components used as-is and are subject to the obligations of their upstream
licenses if the image itself is redistributed.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Notebooks plugin and impose no obligation on users
of the plugin.
