const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  mode: 'development',
  entry: {
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.tsx', '.js', '.json'],
  },
  module: {
    rules: [
      {test: /\.tsx?$/, loader: 'babel-loader', options: {
        'plugins': ['@vue/babel-plugin-jsx'],
      }},
      {test: /\.tsx?$/, loader: 'ts-loader', options: {allowTsInNodeModules: true}},
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader', 'postcss-loader'],
      },
      {
        test: /\.(js|mjs|jsx|ts|tsx)$/,
        enforce: 'pre',
        use: ['source-map-loader'],
      },
    ],
  },
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'DG': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'cash-dom': '$',
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
    'html2canvas': 'html2canvas',
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
    clean: true,
  },
};
