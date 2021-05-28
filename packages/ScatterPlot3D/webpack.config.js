const path = require('path');

module.exports = {
    mode: 'development',
    entry: {
        package: './package.js'
    },
    // devtool: 'inline-source-map',
    externals: {
        'datagrok-api/dg': 'DG',
        'datagrok-api/grok': 'grok',
        'datagrok-api/ui': 'ui',
        "test22": "test22",
        "rxjs": "rxjs",
        "rxjs/operators": "rxjs.operators"
    },
    output: {
        filename: '[name].js',
    //    filename: 'ee.js',
   //     library: 'charts',
     //   libraryTarget: 'var',
        path: path.resolve(__dirname, 'dist'),
        //path: './dist'
    },
};