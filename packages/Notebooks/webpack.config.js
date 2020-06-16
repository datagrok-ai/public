const path = require('path');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');

module.exports = {
    mode: 'development',
    entry: {
        package: './src/package.js'
    },
    //devtool: 'inline-source-map',
    externals: {
        'datagrok-api/dg': 'DG',
        'datagrok-api/grok': 'grok',
        'datagrok-api/ui': 'ui',
        "openchemlib/full.js": "OCL",
        "rxjs": "rxjs",
        "rxjs/operators": "rxjs.operators"
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
            { test: /\.css$/, use: ['style-loader', 'css-loader'] },
            { test: /\.html$/, use: 'file-loader' },
            { test: /\.md$/, use: 'raw-loader' },
            { test: /\.js.map$/, use: 'file-loader' },
            {
                // In .css files, svg is loaded as a data URI.
                test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
                issuer: { test: /\.css$/ },
                use: {
                    loader: 'svg-url-loader',
                    options: { encoding: 'none', limit: 10000 }
                }
            },
            {
                // In .ts and .tsx files (both of which compile to .js), svg files
                // must be loaded as a raw string instead of data URIs.
                test: /\.svg(\?v=\d+\.\d+\.\d+)?$/,
                issuer: { test: /\.js$/ },
                use: {
                    loader: 'raw-loader'
                }
            },
            {
                test: /\.(png|jpg|gif|ttf|woff|woff2|eot)(\?v=[0-9]\.[0-9]\.[0-9])?$/,
                use: [{ loader: 'url-loader', options: { limit: 10000 } }]
            }
        ]
    },
    plugins: [
        new CleanWebpackPlugin()
    ],
};
