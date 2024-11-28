const path = require('path');
const {CleanWebpackPlugin} = require('clean-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'production',
  entry: {
    package: './src/package.js',
    jupyterStyles: [
      '@jupyterlab/application/style/index.css',
      '@jupyterlab/codemirror/style/index.css',
      '@jupyterlab/completer/style/index.css',
      '@jupyterlab/documentsearch/style/index.css',
      '@jupyterlab/notebook/style/index.css',
      './css/application-base.css',
      './css/notebooks.css',
      './css/theme-light-extension-index.css',
      './css/ui-components-base.css',
    ]
  },
  // devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    "openchemlib/full.js": "OCL",
    "rxjs": "rxjs",
    "rxjs/operators": "rxjs.operators"
  },
  resolve: {
    fallback: {
      'path': false,
    }
  },
  optimization: {
    minimize: true
  },
  output: {
    filename: '[name].js',
    library: 'notebooks',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  module: {
    rules: [
      {
        test: /\.css$/,
        exclude: /@jupyterlab.*\.css$/,
        include: /notebooks\.css$/,
        use: [
          'style-loader',
          {
            loader: 'css-loader',
            options: {
              import: {
                filter: (url, media, resourcePath) =>
                !((resourcePath.endsWith('ui-components/style/index.css') && url.includes('base.css')) ||
                  (resourcePath.endsWith('application/style/index.css') && url.includes('base.css')) ||
                  (url.includes('blueprint.css')))
              }
            }
          }
        ]
      },
      {
        test: /\.css$/,
        exclude: /notebooks\.css$/,
        use: [
          MiniCssExtractPlugin.loader,
          'css-loader',
        ]
      },
      {test: /\.html$/, use: 'file-loader'},
      {test: /\.md$/, use: 'raw-loader'},
      {test: /\.js.map$/, use: 'file-loader'},
      {
        // In .css files, svg is loaded as a data URI.
        test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
        issuer: /\.css$/,
        use: {
          loader: 'svg-url-loader',
          options: {encoding: 'none', limit: 10000}
        }
      },
      {
        // In .ts and .tsx files (both of which compile to .js), svg files
        // must be loaded as a raw string instead of data URIs.
        test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
        issuer: /\.js$/,
        use: {
          loader: 'raw-loader'
        }
      },
      {
        test: /\.(png|jpg|gif|ttf|woff|woff2|eot)(\?v=[0-9]\.[0-9]\.[0-9])?$/,
        use: [{loader: 'url-loader', options: {limit: 10000}}]
      }
    ]
  },
  plugins: [
    new CleanWebpackPlugin(),
    new MiniCssExtractPlugin({ filename: 'styles/jupyter-styles.css'})
  ]
};
