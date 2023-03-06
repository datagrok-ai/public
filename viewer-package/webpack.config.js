const fs = require('fs');
const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const fileNames = fs.readdirSync('./src')
  .filter((f) => f !== 'package-test.ts' && f !== 'package.ts' && f.endsWith('.ts'))
  .map((f) => `./src/${f}`);

module.exports = {
  mode: 'development',
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts',
    functions: {filename: 'package-functions.js', library: {type: 'var', name:`${packageName}_functions`}, import: fileNames},
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
  },
  resolveLoader: {
    alias: {
      'func-gen-loader': path.join(path.dirname(__dirname), 'tools', 'loaders', 'func-gen-loader.js'),
    },
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' },
      { test: /(?<!package|package-test)\.tsx?$/, use: [
        'ts-loader', 'func-gen-loader'
      ]}
    ],
  },
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
    'wu': 'wu',
    'exceljs': 'ExcelJS',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
