const path = require('path');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.tsx', '.js'],
  },
  module: {
    rules: [
      {test: /\.js$/, enforce: 'pre', use: ['source-map-loader'], exclude: /node_modules/},
      {test: /\.ts(x?)$/, use: {loader: 'ts-loader', options: {allowTsInNodeModules: true}}},
      {test: /\.css$/, use: ['style-loader', 'css-loader']},
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
  },
  output: {
    filename: '[name].js',
    library: 'hittriage',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
