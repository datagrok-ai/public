# Building the vendored Patinae viewer

`wasm/` holds the built web viewer from https://github.com/zmactep/patinae (`web/`):
`patinae-viewer.js` (TypeScript API), `patinae_web.js` (wasm-bindgen glue, renamed from the hashed
upstream file), `patinae_web_bg.wasm`, and `VERSION` (tag, commit, applied patches). The upstream
panel stylesheet is vendored as `css/patinae-panels.css`, wrapped in `.patinae-view { }`.

Prerequisites: Rust stable with the `wasm32-unknown-unknown` target
(`rustup target add wasm32-unknown-unknown`), Node 22. `wasm-pack` is installed by upstream's
`npm install`.

```
npm run build:wasm            # bash wasm-src/build.sh
PATINAE_TAG=v0.4.8 npm run build:wasm
```

The script clones the tag into a temp dir (or `PATINAE_WORK`), applies `patches/*.patch`, runs
`npm run build` in `web/`, copies the artefacts, and rewrites `wasm/VERSION`. Commit the regenerated
files together.

## Patches

- `0001-web-load-pse.patch` — wires the `pse` format into `WebViewer::load_data` using the
  session crate's `load_pse_bytes` + `pse_to_session` (the native `load` command already does this;
  the web bridge only handled `prs`). Drop it once upstream ships the same change.
- `0002-render-unaligned-object-entry.patch` — `SceneStore::object_entry` read an `ObjectEntry`
  out of a `Vec<u8>` with `bytemuck::from_bytes`, which panics
  (`TargetAlignmentGreaterAndInputNotAligned`) inside `render_frame` on Chrome/Windows and stops
  the viewer with a black canvas; upstream's own prebuilt wasm shows the same panic.
  `pod_read_unaligned` fixes it. Drop once upstream merges the equivalent.
