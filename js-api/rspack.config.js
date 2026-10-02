const path = require('path');
const {swcOptions, loaders} = require('@datagrok/build-config');

const ts = (target) => {
  const options = swcOptions({});
  options.jsc.target = target;
  return {test: /\.ts$/, exclude: /node_modules/, loader: 'builtin:swc-loader', options,
    resolve: {extensionAlias: {'.js': ['.ts', '.js']}}};
};

module.exports = (env = {}) => {
  const common = {
    mode: env.only === 'browser' || process.env.NODE_ENV === 'production' ? 'production' : 'development',
    devtool: 'source-map',
    resolve: {extensions: ['.ts', '.js'], alias: {'cash-dom$': require.resolve('cash-dom/dist/cash.esm.js')}},
  };
  const node = {
    ...common,
    target: 'node',
    entry: {DG: './datagrok.ts'},
    module: {rules: [ts('es2020'), {test: /\.css$/i, loader: loaders.null}]},
    output: {filename: 'datagrok.js', library: {type: 'commonjs'}, path: __dirname},
  };
  const browser = {
    ...common,
    entry: {DG: './dg.ts'},
    module: {rules: [
      ts('es2015'),
      {test: /[\\/](typeahead-standalone|libraries[\\/][^\\/]+[\\/]dist)[\\/].*\.m?js$/, loader: 'builtin:swc-loader',
        options: {env: {targets: 'chrome 50'}, jsc: {parser: {syntax: 'ecmascript'}}}},
      {test: /\.css$/i, use: [loaders.style, {loader: loaders.css, options: {url: false}}], type: 'javascript/auto'},
    ]},
    externals: {'openchemlib/full.js': 'OCL', 'rxjs': 'rxjs', 'rxjs/operators': 'rxjs.operators'},
    output: {filename: 'js-api.js', library: {name: '[name]', type: 'var'}, path: path.resolve(__dirname, '../../core/client/xamgle/web/js/api')},
  };
  return env.only ? [{node, browser}[env.only]] : [node, browser];
};
