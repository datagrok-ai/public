const path = require('path');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');

module.exports = {
    mode: 'production',
    entry: {'grok-api': './src/grok-api.js'},
    //  devtool: 'inline-source-map',
    devServer: {
        contentBase: './dist',
        hot: true,
    },
    module: {
        rules: [
            {
                test: /\.m?js$/,
                exclude: /(node_modules|bower_components)/,
                use: {
                    loader: 'babel-loader',
                    options: {
                        "presets": [
                            ["@babel/preset-env", {
                                "targets": { "browsers": ["last 2 chrome versions", "chrome 50"] },
                                "useBuiltIns": "usage"
                            }]
                        ]
                    }
                }
            }
        ]
    },
    plugins: [],
    output: {
        filename: '[name].js',
        library: 'grok',
        libraryTarget: 'var',
        path: path.resolve(__dirname, '../../xamgle/web/js/api'),
    },
};