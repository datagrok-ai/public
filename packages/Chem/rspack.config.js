const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "openchemlib/full": "OCL",
    "NGL": "NGL"
  },
  "wasm": "async",
  "resolve": {
    "fallback": {
      "perf_hooks": false
    }
  }
});
