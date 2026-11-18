const {bundler, rspack} = require('@datagrok/build-config');

module.exports = bundler({
  "plugins": [
    new rspack.CopyRspackPlugin({patterns: [
      {from: 'crux', to: 'crux', globOptions: {ignore: ['**/README.md', '**/engine.toml']}, info: {minimized: true}},
    ]}),
  ],
});
