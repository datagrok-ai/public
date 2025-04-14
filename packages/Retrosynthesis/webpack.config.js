const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const mode = 'production';
const webpack = require('webpack');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: mode,
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts'
  },
  resolve: {
    fallback: { "url": false },
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' },
      { test: /\.css$/, use: ['style-loader', 'css-loader'] },
      { test: /\.(jpe?g|gif|png|svg|sdf)$/, loader: "file-loader" }
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui'
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new webpack.DefinePlugin({
        process: {env: {}}
    })
  ]
};
