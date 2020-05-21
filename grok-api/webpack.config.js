const path = require('path');

module.exports = {
    mode: 'development',
    entry: {'grok': './grok.js', 'ui': './ui.js', 'DG': './dg.js'},
    devtool: 'inline-source-map',
    module: {
        rules: [
            {
                test: /\.m?js$/,
                exclude: /node_modules/,
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
    externals: {"openchemlib/full.js": "OCL"},
    plugins: [],
    optimization: {
        splitChunks: {
            cacheGroups: {
                vendor: {
                    chunks: 'initial',
                    name: 'common',
                    enforce: true
                },
            }
        }
    },
    output: {
        filename: 'grok-api-[name].js',
        library: '[name]',
        libraryTarget: 'var',
        path: path.resolve(__dirname, '../../xamgle/web/js/api'),
    },
};