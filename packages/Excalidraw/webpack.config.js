const path = require('path');
const webpack = require('webpack');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const mode = 'development';

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: mode,
  entry: {
    package: './src/package.ts',
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
  },
  resolve: {
    fallback: { "url": false },
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
    alias: {
      'roughjs/bin/rough': 'roughjs/bin/rough.js',
      'roughjs/bin/math': 'roughjs/bin/math.js',
      'roughjs/bin/generator': 'roughjs/bin/generator.js',
    }
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
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators'
  },
  output: {
    filename: '[name].js',
    library: 'excalidraw',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new webpack.DefinePlugin({
        process: {env: {}}
    })
  ]
};
