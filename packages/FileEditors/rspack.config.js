const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "package-rtf.js": "RTFJS"
  }
});
