const path = require('path');

module.exports = {
    mode: 'development',
    devtool: 'inline-source-map',

    entry: {
        package: './src/package.ts'
    },
    resolve: {
        extensions: [".ts", ".tsx", ".js", ".jsx", ".css"]
    },

    module: {
        rules: [
            {
                test: /\.ts(x?)$/,
                exclude: /node_modules/,
                use: [
                    {
                        loader: "ts-loader"
                    }
                ]
            },
            {
                enforce: 'pre',
                test: /\.js$/,
                loader: "source-map-loader"
            },
            {
                test: /\.css$/i,
                use: ['style-loader', 'css-loader'],
            },
        ]
    },

    externals: {
        'datagrok-api/dg': 'DG',
        'datagrok-api/grok': 'grok',
        'datagrok-api/ui': 'ui',
        "openchemlib/full.js": "OCL",
        "rxjs": "rxjs",
        "rxjs/operators": "rxjs.operators"
    },

    output: {
        filename: '[name].js',
        library: 'sunburst',
        libraryTarget: 'var',
        path: path.resolve(__dirname, 'dist'),
    },
};
