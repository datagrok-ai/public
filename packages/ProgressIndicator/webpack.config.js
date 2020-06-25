const path = require('path');

module.exports = {
    mode: 'development',
    entry: {
        package: './src/package.js'
    },
    devtool: 'inline-source-map',
    externals: {
        'datagrok-api/dg': 'DG',
        'datagrok-api/grok': 'grok',
        'datagrok-api/ui': 'ui',
        "openchemlib/full.js": "OCL",
        "rxjs": "rxjs",
        "rxjs/operators": "rxjs.operators"
    },
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
    output: {
        filename: '[name].js',
        library: 'progressindicator',
        libraryTarget: 'var',
        path: path.resolve(__dirname, 'dist'),
    },
};