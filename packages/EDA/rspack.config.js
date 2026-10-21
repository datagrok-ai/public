const {bundler, rspack} = require('@datagrok/build-config');

module.exports = bundler({
  "emit": ["./wasm/sci_comp_ml_bg.wasm", "./wasm/XGBoostAPI.wasm", "./wasm/SVMAPI.wasm"],
  "plugins": [
    new rspack.CopyRspackPlugin({patterns: [
      {from: 'crux', to: 'crux', globOptions: {ignore: ['**/README.md']}, info: {minimized: true}},
    ]}),
  ],
});
