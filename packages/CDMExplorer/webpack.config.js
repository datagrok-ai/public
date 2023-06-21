const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts'
  },
  resolve: {
    extensions: ['.js', '.json', '.ts'],
  },
  module: {
    rules: [
      { test: /\.ts$/, loader: 'ts-loader' }
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    "openchemlib/full.js": "OCL",
    "rxjs": "rxjs",
    "rxjs/operators": "rxjs.operators",
    'cash-dom': '$',
    'dayjs': 'dayjs',
  },
  output: {
    filename: '[name].js',
    library: 'cdm',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
