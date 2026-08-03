// The Node.js runtime bundle (datagrok.js) alone — the first target of the main
// config. CI test pipelines build the browser bundle via webpack.config.docker.js
// and add this one so `grok test`'s Node pass can load datagrok-api/datagrok.
module.exports = require('./webpack.config.js')[0];
