const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  jsx: 'react',
  resolve: {
    alias: {
      'roughjs/bin/rough': 'roughjs/bin/rough.js',
      'roughjs/bin/math': 'roughjs/bin/math.js',
      'roughjs/bin/generator': 'roughjs/bin/generator.js',
    },
  },
});
