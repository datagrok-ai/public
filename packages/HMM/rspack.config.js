const {bundler} = require('@datagrok/build-config');

// The engine (wasm/hmmer_web.wasm) and the ANARCI data are fetched at run time
// from dist/ by the workers, so they are copied there rather than bundled.
module.exports = bundler({
  emit: ['wasm/hmmer_web.wasm', 'wasm/anarci/*', 'wasm/pfam/*'],
});
