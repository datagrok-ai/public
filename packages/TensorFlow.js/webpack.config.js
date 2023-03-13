const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.js',
    utils: './src/utils.js',
    nn: './src/nn.js',
    'train.worker': './src/train.worker.js'
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs'
  },
  output: {
    filename: '[name].js',
    library: 'tensorflowdev',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  }
};
