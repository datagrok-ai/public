const path = require('path');

module.exports = {
  mode: 'production',
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name:`demo_test`}, import: './src/package-test.ts'},
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.js', '.json', '.ts'],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' }
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
  },
  output: {
    filename: '[name].js',
    library: 'demo',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
