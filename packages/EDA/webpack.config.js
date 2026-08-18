const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, ''); 
const mode = 'production';
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');


module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: mode,
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
    // The wasm binaries are listed so webpack emits them to dist/
    // (via the .wasm file-loader rule below) for the runtime init URLs.
    package: ['./wasm/sci_comp_ml_bg.wasm', './wasm/XGBoostAPI.wasm', './src/package.ts']
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.js', '.json', '.tsx'],
  },
  module: {
    rules: [
      // Emit the wasm-bindgen binary to dist/ with its original name so the
      // loader can fetch `${webRoot}/dist/sci_comp_ml_bg.wasm` (mirrors Chem).
      {
        test: /\.wasm$/i,
        type: 'javascript/auto',
        loader: 'file-loader',
        options: {publicPath: 'dist/', name: '[name].[ext]'},
      },
      {test: /\.js$/, enforce: 'pre', use: ['source-map-loader'], exclude: /node_modules/},
      {test: /\.ts(x?)$/, use: 'ts-loader', exclude: /node_modules/},
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
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
    'cash-dom': '$',
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})
  ]
};
