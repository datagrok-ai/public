const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "openchemlib/full": "OCL",
    "NGL": "NGL"
  },
  "emit": ["./src/RDKit_minimal_1.2.23.wasm", "./src/bbwasm_bg.wasm", "./src/crux/crux_wasm_bg.wasm"],
  "resolve": {
    "fallback": {
      "perf_hooks": false
    }
  }
});
