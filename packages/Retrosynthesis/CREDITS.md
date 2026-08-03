# Retrosynthesis — Third-Party Libraries

The `@datagrok/retrosynthesis` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

The TypeScript bundle has only one third-party runtime JavaScript dependency
(`yaml`, ISC). The heavy lifting — AiZynthFinder retrosynthetic tree search —
runs server-side inside the Docker image (Section 4).

---

## 1. Bundled in the published artifact (`dist/`)

### yaml (2.x)

YAML 1.1/1.2 parser used to read AiZynthFinder config folders' `config.yml`
files when reconciling expansion / stock / filter policy selections.

- Upstream: https://github.com/eemeli/yaml
- License: **ISC**

```
Copyright Eemeli Aro <eemeli@gmail.com>

Permission to use, copy, modify, and/or distribute this software for any purpose
with or without fee is hereby granted, provided that the above copyright notice
and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

None used at runtime by this package's own code beyond the standard
`datagrok-api/*` namespaces and `rxjs` provided by the platform.

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

These images are built and distributed separately from the JavaScript bundle.

### `dockerfiles/aizynthfinder` — AiZynthFinder Celery worker

The image is built from `python:3.10`, installs `aizynthfinder[all]==4.3.2`
into a wheel cache, then ships the slim Python 3.10 image with the engine
binaries copied in. The `download_public_data` AiZynthFinder helper is run
at image-build time to fetch the default policy / stock / model files into
`/app/configs/default/`.

| Component | Source | License |
|-----------|--------|---------|
| Python 3.10 | https://www.python.org/ | PSF License |
| **AiZynthFinder 4.3.2** | https://github.com/MolecularAI/aizynthfinder | **MIT** (AstraZeneca) |
| Celery | https://github.com/celery/celery | BSD-3-Clause |
| filelock | https://github.com/tox-dev/py-filelock | Unlicense (public domain) |
| NumPy 1.26.4 | https://numpy.org/ | BSD-3-Clause |
| ONNX Runtime | https://github.com/microsoft/onnxruntime | MIT |
| PyTorch | https://pytorch.org/ | BSD-3-Clause-style (PyTorch license) |
| RDKit | https://www.rdkit.org/ | BSD-3-Clause |
| patchelf | https://github.com/NixOS/patchelf | GPL-3.0 (build-time only — not redistributed in the runtime image) |
| datagrok-celery-task | (Datagrok-internal) | covered by repo MIT |
| datagrok-api (Python) | (Datagrok-internal) | covered by repo MIT |

The default AiZynthFinder model / stock / template files downloaded by
`download_public_data` are made available by AstraZeneca / MolecularAI under
the AiZynthFinder MIT license. The container additionally fetches Debian
system libraries (`libgomp1`, `libpq5`, `libxrender1`, `libxext6`, `libsm6`,
`libgl1`) from the upstream Debian distribution under their respective
licenses. `patchelf` is invoked once at image-build time to clear the
executable-stack flag on a single ONNX Runtime shared object (per Microsoft
ONNX Runtime issue 24911) and is **not** redistributed inside the runtime
image.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI, the TypeScript / webpack toolchain, and the
`file-loader` webpack plugin (MIT, *Copyright JS Foundation and other
contributors*) used at build time only. The dependency on `@types/react-dom`
provides TypeScript type definitions only and contributes no runtime code to
the bundle. The peer/devDependency on `@datagrok/chem` is MIT — covered by
the repo-wide `LICENSE.md`. These are **not** redistributed as part of the
published Retrosynthesis plugin and impose no obligation on users of the
plugin.
