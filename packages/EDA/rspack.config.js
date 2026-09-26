const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "emit": ["./wasm/sci_comp_ml_bg.wasm", "./wasm/XGBoostAPI.wasm", "./wasm/SVMAPI.wasm"]
});
