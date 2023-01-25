const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  mode: 'development',
  entry: {
    package: ['./src/arrow1_bg.wasm', './src/package.ts'],
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
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
  devtool: 'inline-source-map',
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