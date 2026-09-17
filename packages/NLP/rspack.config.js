const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "resolve": {
    "fallback": {
      "util": false
    }
  }
});
