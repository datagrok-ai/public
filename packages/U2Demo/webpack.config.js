const path = require('path');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts',
  },
  resolve: {
    symlinks: false,
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
    // build against the library itself, not the copy `file:` deps leave in node_modules
    // (that copy goes stale on every u2 change until a reinstall)
    alias: {'@datagrok-libraries/u2': path.resolve(__dirname, '../../libraries/u2')},
  },
  module: {
    rules: [
      // `?raw` inlines a source file as text so the demo can show the source of the build that is
      // running (src/source-panel.ts); everything else compiles as usual
      {resourceQuery: /raw/, type: 'asset/source'},
      {test: /\.tsx?$/, resourceQuery: {not: [/raw/]}, loader: 'ts-loader',
        options: {allowTsInNodeModules: true}},
      {test: /\.css$/i, use: ['style-loader', 'css-loader']},
    ],
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'datagrok-api/u2core': 'DG.U2',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
    'html2canvas': 'html2canvas',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
