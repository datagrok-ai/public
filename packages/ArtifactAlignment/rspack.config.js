const path = require('path');
const fs = require('fs');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');

// Resolve linked-or-npm-installed lib paths to their real on-disk location so
// module.rules[].include matches via the resolved path (not the node_modules symlink).
const realLib = (pkg) => {
  try { return fs.realpathSync(path.resolve(__dirname, 'node_modules', pkg)); }
  catch { return null; }
};
const SOURCE_LIBS = [
  '@datagrok-libraries/compute-utils',
  '@datagrok-libraries/utils',
  '@datagrok-libraries/test',
].map(realLib).filter(Boolean);

module.exports = (env = {}) => ({
  mode: env.NODE_ENV ?? 'development',
  cache: false,
  entry: {
    test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
    mainFields: ['source', 'browser', 'module', 'main'],
  },
  module: {
    rules: [
      {
        test: /\.ts$/,
        include: [path.resolve(__dirname, 'src'), ...SOURCE_LIBS],
        loader: 'builtin:swc-loader',
        options: {
          jsc: {
            target: 'es2015',
            parser: {syntax: 'typescript', decorators: true, tsx: false},
            transform: {legacyDecorator: true, decoratorMetadata: true},
          },
        },
      },
    ],
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
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
  },
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
    clean: true,
  },
});
