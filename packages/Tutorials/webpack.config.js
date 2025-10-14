const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
  
module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    package: './src/package.ts',
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' },
      { test: /\.css$/i, use: ['style-loader', 'css-loader'] },
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
  },
  output: {
    filename: '[name].js',
    library: 'tutorials',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})
  ]
};
