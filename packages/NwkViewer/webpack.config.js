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
  },
  output: {
    filename: '[name].js',
    library: 'nwkviewer',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
};
