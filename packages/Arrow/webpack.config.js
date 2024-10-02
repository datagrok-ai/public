const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  experiments: {
    syncWebAssembly: true,
  },
  entry: {
    package: './src/package.ts',
    test: {
      filename: 'package-test.js',
      library: { type: 'var', name: `${packageName}_test` },
      import: './src/package-test.ts',
    },
  },
  resolve: {
    extensions: ['.wasm', '.ts', '.mjs', '.js', '.json', '.tsx'],
  },
  module: {
    rules: [
      {
        test: /\.(wasm)$/i,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/",
          name: '[name].[ext]'
        }
      },
      {
        test: /\.ts?$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
  },
  output: {
    filename: '[name].js',
    library: 'arrow',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
