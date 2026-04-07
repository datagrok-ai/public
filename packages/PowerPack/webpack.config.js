const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
const workingInDartium = true; // when working in Dartium, set this to true

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    package: './src/package.ts',
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'}
  },
  resolve: {
    extensions: ['.ts', '.wasm', '.mjs', '.json', '.js', '.tsx'],
  },
  module: {
    rules: [
      { test: /\.ts$/,
        loader: 'ts-loader',
        exclude: /node_modules/,
      },
      ...(workingInDartium ? [{
        test: /\.js$/,
        // include: [
        //   /node_modules[\\/]konva/,
        //   /node_modules[\\/]@datagrok-libraries/,
        // ],
        use: {
          loader: 'babel-loader',
          options: {
            presets: [['@babel/preset-env', {
              targets: 'chrome 54'// approximate Dartium's engine
            }]],
          },
        },
      }] : []),
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
  },
  output: {
    filename: '[name].js',
    library: 'powerpack',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
};
