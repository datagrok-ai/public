const path = require('path');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

const mode = process.env.NODE_ENV ?? 'production';
if (mode !== 'production')
  console.warn(`Building 'BiostructureViewer' in '${mode}' mode.`);

module.exports = {
  ...(mode === 'production' ? {cache: {type: 'filesystem'}} : {}),
  mode: mode,
  entry: {
    package: ['./src/package.ts'],
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
  },
  resolve: {
    extensions: ['.ts', '.tsx', '.wasm', '.mjs', '.js', '.json'],
    fallback: {
      'url': false,
      'fs': false,
      'path': false, //require.resolve('path-browserify'),
      'crypto': false, // require.resolve('crypto-browserify'),
    },
  },
  module: {
    rules: [
      {test: /\.js$/, enforce: 'pre', use: ['source-map-loader'], exclude: /node_modules/},
      {test: /\.ts(x?)$/, loader: 'ts-loader', exclude: /node_modules/},
      {
        test: /\.(html|ico)$/,
        use: [
          {loader: 'file-loader', options: {name: '[name].[ext]'}},
        ],
      },
      {
        test: /\.(s*)css$/,
        use: [
          'style-loader',
          MiniCssExtractPlugin.loader,
          {loader: 'css-loader', options: {sourceMap: false, modules: {localIdentName: '[local]'}}},
          {loader: 'sass-loader', options: {sourceMap: false}},
        ],
      },
    ],
  },
  plugins: [
    new MiniCssExtractPlugin({filename: 'molstar.css'}),
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
  devtool: mode !== 'production' ? 'inline-source-map' : 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
    'NGL': 'NGL',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
