const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.js'
  },
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
                "targets": {"browsers": ["last 2 chrome versions", "chrome 50"]},
                "useBuiltIns": "usage"
              }]
            ]
          }
        }
      }
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
    library: 'impute',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};