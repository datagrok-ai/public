const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "externals": {
    "dayjs/plugin/utc": "utc"
  }
});
