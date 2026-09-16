// Compatibility marker for stands that predate the workspace toolchain: the server marks a
// published version as bundled (isWebpack) only when webpack.config.js is present in the
// archive, and the client requests dist/package.js only for bundled versions. The real
// build is `grok build` (rspack.config.js); devtool and externals below only satisfy
// `grok check` and mirror the rspack config.
module.exports = {
  devtool: 'source-map',
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
    'wu': 'wu',
    'exceljs': 'ExcelJS',
    'html2canvas': 'html2canvas',
    'vue': 'Vue',
  },
};
