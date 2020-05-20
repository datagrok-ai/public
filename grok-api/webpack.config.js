const path = require('path');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');

module.exports = {
    mode: 'production',
    entry: './src/grok-api.js',
  //  devtool: 'inline-source-map',
    devServer: {
        contentBase: './dist',
        hot: true,
    },
    plugins: [],
    output: {
        filename: 'grok-api.js',
        library: 'grok',
        libraryTarget: 'var',
        path: path.resolve(__dirname, '../../xamgle/web/js/api'),
    },
};