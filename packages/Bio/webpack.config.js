const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts',
  },
  resolve: {
    fallback: {'url': false},
    extensions: ['.wasm', '.mjs', '.ts', '.js', '.json', '.tsx'],
  },
  module: {
    rules: [
      {
        test: /\.ts(x?)$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader'],
        exclude: /node_modules/,
      },
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
    library: 'bio',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
