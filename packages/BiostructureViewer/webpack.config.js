const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts',
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
    fallback: {
      'fs': false,
      'path': false, //require.resolve('path-browserify'),
      'crypto': false, // require.resolve('crypto-browserify'),
    },
  },
  module: {
    rules: [
      {test: /\.tsx?$/, loader: 'ts-loader'},
      {
        test: /\.(html|ico)$/,
        use: [{
          loader: 'file-loader',
          options: {name: '[name].[ext]'},
        }],
      },
      {
        test: /\.(s*)css$/,
        use: [
          MiniCssExtractPlugin.loader,
          {loader: 'css-loader', options: {sourceMap: false}},
          {loader: 'sass-loader', options: {sourceMap: false}},
        ],
      },
    ],
  },
  plugins: [
    new MiniCssExtractPlugin({filename: 'molstar.css'}),
  ],
  devtool: 'inline-source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
  },
  output: {
    filename: '[name].js',
    library: 'biostructureviewer',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
