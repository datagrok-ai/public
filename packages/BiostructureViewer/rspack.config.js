const {bundler, rspack, loaders} = require('@datagrok/build-config');

// Mol* styles are extracted to molstar.css; class names are kept verbatim (localIdentName: [local]).
module.exports = bundler({
  jsx: 'react',
  externals: {NGL: 'NGL'},
  css: false,
  rules: [
    {test: /\.(html|ico)$/, type: 'asset/resource', generator: {filename: '[name][ext]'}},
    {
      test: /\.(s*)css$/,
      use: [
        rspack.CssExtractRspackPlugin.loader,
        {loader: loaders.css, options: {sourceMap: false, modules: {localIdentName: '[local]'}}},
        {loader: 'sass-loader', options: {sourceMap: false}},
      ],
      type: 'javascript/auto',
    },
  ],
  plugins: [new rspack.CssExtractRspackPlugin({filename: 'molstar.css'})],
});
