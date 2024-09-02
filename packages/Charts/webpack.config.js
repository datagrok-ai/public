const path = require('path');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts'
  },
  resolve: {
    extensions: ['.ts', '.json', '.js'],
  },
  module: {
    rules: [
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
      {
        test: /\.ts$/,
        loader: 'ts-loader',
      },
    ]
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'wu': 'wu',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
