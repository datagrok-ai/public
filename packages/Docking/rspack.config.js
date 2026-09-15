const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "NGL": "NGL"
  }
});
