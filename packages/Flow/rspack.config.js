const {bundler} = require('@datagrok/build-config');

// rete-history-plugin optionally imports rete-comment-plugin, which Flow does not use.
module.exports = bundler({jsx: 'react', resolve: {fallback: {'rete-comment-plugin': false}}});
