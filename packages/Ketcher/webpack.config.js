const path = require('path');
const webpack = require('webpack');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts'
  },
  resolve: {
    fallback: { "url": false },
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
  },

  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' },
      { test: /\.css$/, use: ['style-loader', 'css-loader'] },
      { test: /\.(jpe?g|gif|png|svg|sdf)$/, loader: "file-loader" }
    ],
  },
  devtool: 'inline-source-map',
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
    library: 'ketcher',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new webpack.DefinePlugin({
        process: {env: {}}
    })
  ]
};
