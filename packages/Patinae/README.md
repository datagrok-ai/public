# Patinae

`Patinae` embeds the [Patinae](https://github.com/zmactep/patinae) molecular viewer (Rust, BSD-3,
compiled to WebAssembly, rendered with WebGPU) into [Datagrok](https://datagrok.ai).

## File previews

Selecting a file with one of these extensions in the file browser shows it in Patinae, with its
command line (top), object list (right) and sequence panel (bottom):

| Extension | Content                                                                |
|-----------|------------------------------------------------------------------------|
| `.pse`    | PyMOL session (objects, representations, colors, camera)               |
| `.prs`    | Patinae native session                                                 |
| `.pml`    | PyMOL script; `load` and `@` lines resolve against the script's folder |

Commands typed into the panel use PyMOL syntax (`show cartoon`, `color marine, chain A`,
`select site, byres around 5 ligand`, ...).

Requirements: a browser with WebGPU and CSS nesting (Chrome or Edge 120+). The wasm module (about 4 MB) is
fetched once per session on the first preview.

Limitations: `fetch <id>` and `load https://...` inside scripts or the command line request
RCSB directly from the browser; sessions cannot be saved back; the platform's 10 MB preview
limit applies to the previewed file.

## Rebuilding the viewer

The vendored viewer in `wasm/` and its panel stylesheet `css/patinae-panels.css` are built from
the upstream tag recorded in `wasm/VERSION` plus the patches in `wasm-src/patches/`. See
`wasm-src/README.md`.

## Samples (`files/samples/`, deployed to `System:AppData/Patinae/samples`)

| File                          | What it shows                                                                                                   |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------|
| `hiv-protease-tour.pml`       | 1HSG with indinavir: sticks, labels, distance measurements, three named scenes, a 180-frame keyframed camera movie |
| `trp-cage-ensemble.pml`       | 1L2Y NMR ensemble: 38 states played as a movie while the view rocks                                              |
| `crambin-density.pml`         | 1CRN with a density map contoured as a mesh (`isomesh`) and crystallographic symmetry mates (`symexp`)           |
| `protease-superposition.pml`  | 1HSG and 1OHR superposed with `align`, inhibitors compared in the shared pocket                                  |
| `hiv-protease.pse`            | PyMOL 3 session of 1HSG (representations, colors, selections, stored views)                                     |
| `trp-cage-ensemble.pse`       | PyMOL 3 session of the 1L2Y ensemble (38 states)                                                                |
| `1crn.pse`, `1crn.prs`, `demo.pml` | Minimal crambin sessions and script                                                                        |

`crambin-density.ccp4` is a Gaussian pseudo-density computed by PyMOL's `map_new`, so the sample
needs no external map download. Sessions were written by PyMOL 3.x open-source; Patinae's `.pse`
importer brings molecules, states, colors, selections, settings and stored views, but not movies,
scenes or measurement objects, which is why the animated material lives in the scripts.
Playback note: camera-only movies run at the renderer frame rate; frames that switch coordinate
states cost roughly 200 ms each in Patinae, so the ensemble plays at about 5 fps.
