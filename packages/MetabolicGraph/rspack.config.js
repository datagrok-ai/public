const {bundler, rspack} = require('@datagrok/build-config');

// glpk-wasm loads its .wasm by URL next to the bundle; escher is a page-level global.
module.exports = bundler({
  externals: {escher: 'escher'},
  rules: [
    {test: /escher\.js$/, type: 'asset/resource', generator: {emit: false}},
  ],
  plugins: [
    new rspack.CopyRspackPlugin({patterns: [{from: 'node_modules/glpk-wasm/dist/glpk.all.wasm', to: '[name][ext]'}]}),
  ],
});
