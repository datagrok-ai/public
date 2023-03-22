const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts'
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.tsx'],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' }
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators'
  },
  output: {
    filename: '[name].js',
    library: 'hittriage',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
