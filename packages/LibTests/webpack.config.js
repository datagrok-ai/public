// Compatibility marker for stands that predate the workspace toolchain: the server marks a
// published version as bundled (isWebpack) only when webpack.config.js is present in the
// archive, and the client requests dist/package.js only for bundled versions. The real
// build is `grok build` (rspack via @datagrok/build-config); devtool below only satisfies
// `grok check`'s source-map regex.
module.exports = {devtool: 'source-map'};
