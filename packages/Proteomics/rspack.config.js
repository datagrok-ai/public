const {bundler} = require('@datagrok/build-config');

// tools/*.sh and *.sql are embedded as strings (see src/global.d.ts).
module.exports = bundler({rules: [{test: /\.(sh|sql)$/, type: 'asset/source'}]});
