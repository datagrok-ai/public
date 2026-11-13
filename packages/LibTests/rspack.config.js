const path = require('path');
const {bundler} = require('@datagrok/build-config');

// rxjs is bundled (not external) and pinned to this package's copy: TestScheduler must share one
// rxjs instance with the code under test. Development mode keeps rxjs's /*@__PURE__*/ mixin calls
// that a minifier would drop. No chunk splitting so the fitting worker stays self-contained.
module.exports = bundler({
  mode: 'development',
  externals: {'rxjs': false, 'rxjs/operators': false, 'DG': 'DG'},
  resolve: {alias: {
    'rxjs': path.resolve(__dirname, 'node_modules/rxjs'),
  }},
  optimization: {concatenateModules: false, splitChunks: false, runtimeChunk: false, chunkIds: 'deterministic'},
  // webcomponents-vue emits its .tsx components as Vue JSX; same transform as Compute2
  rules: [{test: /\.jsx$/, loader: 'babel-loader', options: {
    babelrc: false,
    configFile: false,
    plugins: [['@vue/babel-plugin-jsx', {isCustomElement: (tag) => tag.startsWith('dg-') || tag === 'dock-spawn-ts'}]],
  }}],
});
