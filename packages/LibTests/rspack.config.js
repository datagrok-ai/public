const path = require('path');
const {bundler} = require('@datagrok/build-config');

// rxjs is bundled (not external) and pinned to this package's copy: TestScheduler must share one
// rxjs instance with the code under test. Development mode keeps rxjs's /*@__PURE__*/ mixin calls
// that a minifier would drop. No chunk splitting so the fitting worker stays self-contained.
// vue is pinned for the same reason: tests and webcomponents-vue code under test must share one runtime.
module.exports = bundler({
  mode: 'development',
  externals: {'rxjs': false, 'rxjs/operators': false, 'vue': false, 'DG': 'DG'},
  resolve: {alias: {
    'rxjs': path.resolve(__dirname, 'node_modules/rxjs'),
    'vue': path.resolve(__dirname, 'node_modules/vue'),
  }},
  optimization: {concatenateModules: false, splitChunks: false, runtimeChunk: false, chunkIds: 'deterministic'},
});
