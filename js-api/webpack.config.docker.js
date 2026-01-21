const path = require('path');

module.exports = [
 {
    mode: 'production',
    entry: {'DG': './dg.ts'},
    devtool: 'source-map',
    resolve: {
      extensions: ['.ts', '.tsx', '.js']
    },
    module: {
      rules: [
        {
          test: /\.js$/,
          enforce: 'pre',
          use: ['source-map-loader'],
        },
        {
          test: /\.ts(x?)$/,
          exclude: /node_modules/,
          use: [
            {
              loader: 'ts-loader'
            }
          ]
        },
        {
          test: /\.css$/i,
          use: [{
            loader: 'style-loader',
          }, {
            loader: 'css-loader',
            options: {
              url: false
            }
          },]
        }
      ]
    },
    externals: {'openchemlib/full.js': 'OCL', 'rxjs': 'rxjs', 'rxjs/operators': 'rxjs.operators'},
    plugins: [],
    output: {
      filename: 'js-api.js',
      library: '[name]',
      libraryTarget: 'var',
      path: path.resolve(__dirname, '../../core/client/xamgle/web/js/api'),
    },
}];
