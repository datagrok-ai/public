const {bundler} = require('@datagrok/build-config');

// CodeMirror 6 is bundled; the platform serves CodeMirror 5.
module.exports = bundler({
  externals: {
    'codemirror': false,
  },
});
