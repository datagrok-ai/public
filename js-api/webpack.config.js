const path = require('path');

module.exports = {
  mode: 'development',
  entry: {'grok': './grok.ts', 'ui': './ui.ts', 'DG': './dg.ts'},
  devtool: 'inline-source-map',
  resolve: {
    extensions: ['.tsx', '.ts', '.js'],
  },
  module: {
    rules: [
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
      /*      {
        test: /\.m?js$/,
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
          options: {
            'presets': [
              ['@babel/preset-env', {
                'targets': {'browsers': ['last 2 chrome versions', 'chrome 50']},
                'useBuiltIns': 'usage'
              }]
            ]
          }
        }
      }*/
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