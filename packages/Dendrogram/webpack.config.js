const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');


module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'production',
  entry: {
    package: ['./src/package.ts'],
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
  },
  devServer: {
    contentBase: './dist',
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
  resolve: {
    fallback: {'url': false},
    extensions: ['.wasm', '.mjs', '.ts', '.js', '.json', '.tsx'],
  },
  module: {
    rules: [
      {
        test: /\.worker\.ts$/,
        loader: 'worker-loader',
        options: {
          inline: 'fallback', // this creates a separate file
        },
      },
      {
        test: /\.(wasm)$/i,
        type: 'javascript/auto',
        loader: 'file-loader',
        options: {
          publicPath: 'dist/',
          name: '[name].[ext]',
        },
      },
      {
        test: /\.ts(x?)$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader'],
        exclude: /node_modules/,
      },
      {
        test: /\.js$/,
        enforce: 'pre',
        use: ['source-map-loader'],
        exclude: /node_modules/,
      },
    ],
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
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})
  ]
};
