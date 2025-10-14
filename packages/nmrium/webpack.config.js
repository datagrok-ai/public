const path = require('path');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx', '.jsx'],
    alias: {
      './lib/components/StructureEditor': './lib/components/StructureEditor.js',
      './lib/components/SmilesSvgRenderer': './lib/components/SmilesSvgRenderer.js',
      './lib/components/MolfileSvgRenderer': './lib/components/MolfileSvgRenderer.js',
      './lib/components/IdcodeSvgRenderer': './lib/components/IdcodeSvgRenderer.js',
      'cheminfo-font/lib-react-cjs/lib-react-tsx/nmr/Peaks': path.resolve(__dirname, 'node_modules/cheminfo-font/lib-react-cjs/lib-react-tsx/nmr/Peaks.js'),

    },
  },
  module: {
    rules: [
      {test: /\.tsx?$/, loader: 'ts-loader'},
      {test: /\.css$/, use: ['style-loader', 'css-loader']},
    ],
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
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
  },
};
