const path = require('path');
const CopyWebpackPlugin = require('copy-webpack-plugin');

const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: false,
  mode: 'production',
  target: 'web',
  entry: {
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
    // Prevent bundling Node built-ins for browser; some deps (e.g., glpk-wasm) reference 'fs' conditionally
    fallback: {
      fs: false,
    },
  },
  module: {
    rules: [
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
        test: /\.tsx?$/,
        loader: 'ts-loader',
        include: [
          path.resolve(__dirname, 'src'),
        ],
        options: {
          transpileOnly: true,
          compilerOptions: {
            verbatimModuleSyntax: true,
            isolatedModules: true,
          },
        },
      },
      {
        test: /escher\.js$/,
        type: 'asset/resource',
        generator: {
          emit: false,
        },
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
    'html2canvas': 'html2canvas',
    'escher': 'escher',
  },
  plugins: [
    // Ensure glpk-wasm runtime .wasm files are copied next to bundles
    new CopyWebpackPlugin({
      patterns: [
        {
          from: 'node_modules/glpk-wasm/dist/glpk.all.wasm',
          to: '[name][ext]'
        },
      ],
    }),
  ],
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
};
