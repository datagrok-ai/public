const path = require('path');

module.exports = [
  {
    target: 'node',
    mode: 'development',
    entry: {'DG': './datagrok.ts'},
    devtool: 'source-map',
    resolve: {
      extensions: ['.ts', '.tsx', '.js']
    },
    module: {
      rules: [
        {
          test: /\.m?js$/,
          exclude: {
            and: [/node_modules/],
            not: [
              /typeahead/
            ]
          },
          use: {
            loader: 'babel-loader',
            options: {
              'presets': [
                ['@babel/preset-env', {
                  targets: {
                    node: 'current',
                  },
                }]
              ]
            }
          }
        },
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
          use: 'null-loader'
        }
      ]
    },
    //externals: {'cash-dom': '$'},
    plugins: [],
    output: {
      filename: 'datagrok.js',
      libraryTarget: 'commonjs',
      path: path.resolve(__dirname, ''),
    },
  },{
    mode: 'development',
    entry: {'DG': './dg.ts'},
    devtool: 'source-map',
    resolve: {
      extensions: ['.ts', '.tsx', '.js']
    },
    module: {
      rules: [
        {
          test: /\.m?js$/,
          exclude: {
            and: [/node_modules/],
            not: [
              /typeahead/
            ]
          },
          use: {
            loader: 'babel-loader',
            options: {
              'presets': [
                ['@babel/preset-env', {
                  'targets': {'browsers': ['last 2 chrome versions', 'chrome 50']}
                }]
              ]
            }
          }
        },
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
