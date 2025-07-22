const {merge} = require('webpack-merge')
const common = require('./webpack.common.js')

module.exports = merge(common, {
  mode: 'development',
  entry: './dev-server/index.js',
  devtool: 'source-map',
  output: {
    filename: 'bundle.js'
  },
  devServer: {
    static: [
      {
        directory: './dev-server'
      },
      {
        directory: '.'  // Serve files from project root
      }
    ],
    open: true,
    port: 7621
  }
})
