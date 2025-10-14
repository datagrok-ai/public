
const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');
  
const mode = process.env.NODE_ENV ?? 'production';
if (mode !== 'production')
  console.warn(`Building '${packageName}' in '${mode}' mode.`);

module.exports = {
  // ...(mode !== 'production' ? {} : {cache: {type: 'filesystem'}}),
  cache: {type: 'filesystem'},
  mode: mode,
  entry: {
    package: './src/package.ts',
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
    dojo: {
      filename: 'package-dojo.js',
      library: {type: 'var', name: `${packageName}_dojo`},
      import: './helm/dojo/package.ts',
    },
  },
  resolve: {
    fallback: {'url': false},
    extensions: ['.ts', '.tsx', '.js', '.wasm', '.mjs', '.json'],
    alias: {
      'vendor/helm-web-editor': mode === 'production' ?
        path.resolve(__dirname, 'vendor', 'helm-web-editor.production.js') :
        path.resolve(__dirname, 'vendor', 'helm-web-editor.development.js'),
    },
  },
  devServer: {
    contentBase: './dist',
  },
  // amd: {toUrlUndefined: true},
  module: {
    noParse: /vendor/,
    rules: [
      {test: /\.js$/, enforce: 'pre', use: ['source-map-loader'], exclude: [/node_modules/]},
      {test: /\.ts(x?)$/, use: 'ts-loader', exclude: [/node_modules/]},
      {test: /\.css$/, use: ['style-loader', 'css-loader'], exclude: [/node_modules/]},
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
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})
  ]
};
