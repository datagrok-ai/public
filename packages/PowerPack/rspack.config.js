const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "codemirror": false,
    "exceljs": false
  }
});
