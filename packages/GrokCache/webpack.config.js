const path = require('path');
const webpack = require('webpack');

module.exports = {
  mode: 'development',
  entry: {
    package: ['./src/grokcache.wasm', './src/package.js']
  },
  devtool: 'inline-source-map',
  devServer: {
    contentBase: './dist'
  },

  module: {
    rules: [
      // WASM files should not be processed but just be emitted and we want
      // to have their public URL.
      // file-loader -> options -> name =  '[name].[ext]'  - prevents webpack name change
      {
        test: /\.(wasm)$/i,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/",
          name: '[name].[ext]'
        }
      }
    ]
  },
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
  },
  output: {
    filename: '[name].js',
    library: 'grokcache',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  }

};
