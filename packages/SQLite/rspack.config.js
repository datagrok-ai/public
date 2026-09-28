const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  "emit": ["./src/sql-wasm.wasm"]
});
