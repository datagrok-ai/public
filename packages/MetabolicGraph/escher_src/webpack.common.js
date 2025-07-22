const webpack = require('webpack')
const package = require('./package.json')

module.exports = {
  resolve: {
    extensions: ['.ts', '.tsx', '.js', '.jsx', '.json']
  },
  devtool: 'source-map',
  module: {
    rules: [
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader']
      },
      {
        test: /\.(ts|tsx)$/,
        exclude: /node_modules/,
        use: {
          loader: 'ts-loader',
          options: {
            transpileOnly: true,
            compilerOptions: {
              jsx: 'react',
              module: 'esnext',
              moduleResolution: 'node'
            }
          }
        }
      },
      {
        test: /\.jsx?$/,
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
          options: {
            presets: [
              ['@babel/preset-env', {
                targets: {
                  browsers: ['last 2 versions', 'ie >= 11']
                },
                modules: 'auto',
                useBuiltIns: 'usage',
                corejs: 3
              }],
              '@babel/preset-react'
            ],
            plugins: [
              ['@babel/plugin-transform-runtime', {
                corejs: 3,
                helpers: true,
                regenerator: true
              }]
            ]
          }
        }
      },
      // Embed font Definitions
      {
        test: /\.svg$/,
        include: /icons\/font/,
        type: 'asset',
        parser: {
          dataUrlCondition: {
            maxSize: 65000
          }
        },
        generator: {
          filename: 'fonts/[name][ext]'
        }
      },
      {
        test: /\.woff$/,
        include: /icons\/font/,
        type: 'asset',
        parser: {
          dataUrlCondition: {
            maxSize: 65000
          }
        },
        generator: {
          filename: 'fonts/[name][ext]',
          mimetype: 'application/font-woff'
        }
      },
      {
        test: /\.woff2$/,
        include: /icons\/font/,
        type: 'asset',
        parser: {
          dataUrlCondition: {
            maxSize: 65000
          }
        },
        generator: {
          filename: 'fonts/[name][ext]',
          mimetype: 'application/font-woff2'
        }
      },
      {
        test: /\.[ot]tf$/,
        include: /icons\/font/,
        type: 'asset',
        parser: {
          dataUrlCondition: {
            maxSize: 65000
          }
        },
        generator: {
          filename: 'fonts/[name][ext]',
          mimetype: 'application/octet-stream'
        }
      },
      {
        test: /\.eot$/,
        include: /icons\/font/,
        type: 'asset',
        parser: {
          dataUrlCondition: {
            maxSize: 65000
          }
        },
        generator: {
          filename: 'fonts/[name][ext]',
          mimetype: 'application/vnd.ms-fontobject'
        }
      }
    ]
  },
  plugins: [
    new webpack.DefinePlugin({
      ESCHER_VERSION: JSON.stringify(package.version)
    })
  ]
}
