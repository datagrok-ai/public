const path = require('path');
const {bundler, rspack, loaders} = require('@datagrok/build-config');

// Jupyter styles are extracted to styles/jupyter-styles.css (served into the notebook iframe);
// the package's own notebooks.css is injected like everywhere else.
module.exports = bundler({
  mode: 'production',
  css: false,
  resolve: {
    alias: {
      // marked v4 (security override) has no default export; @jupyterlab/rendermime@2.x expects one.
      'marked$': path.resolve(__dirname, 'src/marked-shim.js'),
      'marked-esm$': path.join(path.dirname(require.resolve('marked/package.json')), 'lib', 'marked.esm.js'),
    },
  },
  rules: [
    {test: /notebooks\.css$/, use: [loaders.style, loaders.css], type: 'javascript/auto'},
    {test: /\.css$/, exclude: /notebooks\.css$/, use: [rspack.CssExtractRspackPlugin.loader, loaders.css], type: 'javascript/auto'},
    {test: /\.html$/, type: 'asset/resource'},
    {test: /\.md$/, type: 'asset/source'},
    {test: /\.js\.map$/, type: 'asset/resource'},
    {test: /\.svg(\?v=\d+\.\d+\.\d+)?$/, issuer: /\.css$/, type: 'asset/inline'},
    {test: /\.svg(\?v=\d+\.\d+\.\d+)?$/, issuer: /\.[jt]sx?$/, type: 'asset/source'},
    {test: /\.(png|jpg|gif|ttf|woff|woff2|eot)(\?v=[0-9]\.[0-9]\.[0-9])?$/, type: 'asset', parser: {dataUrlCondition: {maxSize: 10000}}},
  ],
  plugins: [
    new rspack.CssExtractRspackPlugin({filename: 'styles/jupyter-styles.css'}),
    new rspack.DefinePlugin({'process.argv': '[]'}),
  ],
});
