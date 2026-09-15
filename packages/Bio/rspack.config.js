const path = require('path');
const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  externals: {'openchemlib/full': 'OCL'},
  wasm: 'async',
  resolve: {alias: {'./immunum_bg.js': path.resolve(__dirname, 'node_modules/immunum/immunum.js')}},
});
