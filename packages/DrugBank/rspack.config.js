const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "openchemlib/full": "OCL"
  }
});
