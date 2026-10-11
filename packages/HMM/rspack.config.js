const {bundler} = require('@datagrok/build-config');

// The engine (wasm/hmmer_web.wasm), the tools module (wasm/hmmer_tools.wasm:
// alignment, query models, logos; fetched on first use) and the ANARCI data are fetched at run time
// from dist/ by the workers, so they are copied there rather than bundled.
module.exports = bundler({
  emit: ['wasm/hmmer_web.wasm', 'wasm/hmmer_tools.wasm', 'wasm/anarci/*', 'wasm/pfam/*'],
});
