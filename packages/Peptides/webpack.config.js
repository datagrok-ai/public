const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: ['./src/package.ts'],
  },
  devServer: {
    contentBase: './dist',
  },
  target: 'web',
  module: {
    rules: [
      {
        test: /\.ts?$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
    ],
  },
  resolve: {
    extensions: ['.mjs', '.js', '.json', '.ts', '.tsx'],
  },
  devtool: 'inline-source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
  },
  output: {
    filename: '[name].js',
    library: 'peptides',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
};
