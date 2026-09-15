const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "datagrok-api/u2core": "DG.U2"
  }
});
