const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  emit: [
    './vendor/crux/crux_wasm_bg.wasm',
    './vendor/crux/inchi_wasm_bg.wasm',
    './vendor/crux/NOTICE.crux-wasm',
    './vendor/crux/NOTICE.inchi-wasm',
    './vendor/crux/VERSION.crux',
  ],
});
