const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "openchemlib/full": "OCL",
    "NGL": "NGL"
  },
  "emit": [
    "./src/RDKit_minimal_1.2.23.wasm",
    "./src/bbwasm_bg.wasm",
    "./vendor/chem-wasm/chem_wasm_bg.wasm",
    "./vendor/chem-wasm/NOTICE.chem-wasm",
    "./vendor/chem-wasm/VERSION.chem-wasm",
    "./vendor/crux-sketch/NOTICE.sketch-wasm",
    "./vendor/crux-sketch/NOTICE.cdxml-wasm",
    "./vendor/crux-sketch/NOTICE.inchi-wasm"
  ],
  "resolve": {
    "fallback": {
      "perf_hooks": false
    }
  }
});
