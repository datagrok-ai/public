const path = require('path');

module.exports = {
  mode: 'development',
  entry: {'grok': './grok.ts', 'ui': './ui.ts', 'DG': './dg.ts'},
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
  optimization: {
    splitChunks: {
      cacheGroups: {
        api: {
          chunks: 'initial',
          name: 'common',
          enforce: true
        }
      }
    }
  },
  output: {
    filename: 'js-api-[name].js',
    library: '[name]',
    libraryTarget: 'var',
    path: path.resolve(__dirname, '../../xamgle/web/js/api'),
  },
};
