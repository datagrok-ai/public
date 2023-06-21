const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.js'
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators'
  },
  output: {
    filename: '[name].js',
    library: 'multiplot',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
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
      }
    ]
  },

  optimization: {
    moduleIds: 'deterministic',
    chunkIds: 'deterministic'
  },
  resolve: {
    extensions: ['.ts', '.js', '.json']
  }
};
