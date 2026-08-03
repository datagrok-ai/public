# PubChemApi — Third-Party Libraries

The `@datagrok/pubchem-api` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

This package is a thin REST client for the PubChem PUG API. It has no
third-party runtime JavaScript dependencies of its own beyond the
platform-provided webpack externals.

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either Datagrok-internal
(`@datagrok-libraries/*`, `datagrok-api`) and covered by the repo-wide MIT
`LICENSE.md`, or webpack externals provided by the platform host (Section 2).

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

None used at runtime by this package's own code.

---

## 3. Fetched at runtime from third-party data services (not bundled, not redistributed)

### PubChem PUG REST API

This plugin queries the public **PubChem** PUG REST endpoints
(`https://pubchem.ncbi.nlm.nih.gov/rest/pug` and `pug_view`) at runtime via
`grok.dapi.fetchProxy` to perform substructure / similarity / identity
search and to render compound information panels. The PubChem service is
operated by the U.S. National Center for Biotechnology Information (NCBI) and
governed by NCBI's data-use policies; no PubChem data, code, or assets are
bundled into the plugin or redistributed by the Datagrok platform.

---

## 4. Docker container images (`dockerfiles/`)

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack toolchain. The
peer/devDependency on `@datagrok/chem` is MIT — covered by the repo-wide
`LICENSE.md`. These are **not** redistributed as part of the published
PubChemApi plugin and impose no obligation on users of the plugin.
