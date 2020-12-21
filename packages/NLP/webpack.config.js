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
  resolve: {
    fallback: {
      util: require.resolve("util/")
    }
  },
  module: {
    rules: [
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
    ]
  },
  output: {
    filename: '[name].js',
    library: 'nlp',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};