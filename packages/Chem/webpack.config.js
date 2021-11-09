const path = require('path');

module.exports = (env, options) => ({
  mode: 'development',
  entry: {
    package: ['./src/RDKit_minimal_2021.03_17.wasm', './src/package.ts']
  },
  devtool: options.mode !== 'production' ? 'inline-source-map' : 'source-map',
  devServer: {
    contentBase: './dist'
  },
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
    library: 'chem',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true
  },
  resolve: {
    fallback: { "url": false },
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx']
  },
  module: {
    rules: [
      // prevent a webpack name change
      {
        test: /\.(wasm)$/i,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/",
          name: '[name].[ext]'
        }
      },
      {
        test: /.worker\.js$/i,
        loader: "worker-loader",
        options: {
          inline: "no-fallback"
        },
      },
      {
        test: /\.tsx?$/,
        loader: 'ts-loader'
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader']
      }
    ],
  }
});