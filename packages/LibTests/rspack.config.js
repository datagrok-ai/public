const path = require('path');
const fs = require('fs');
const rspack = require('@rspack/core');
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
  '@datagrok-libraries/compute-api',
  '@datagrok-libraries/utils',
  '@datagrok-libraries/webcomponents',
  '@datagrok-libraries/test',
  '@datagrok-libraries/arrow',
  'diff-grok',
].map(realLib).filter(Boolean);

module.exports = (env = {}) => {
  // Development mode (no minification) — minifiers drop calls marked /*@__PURE__*/
  // by rxjs's _esm5 build (e.g. `applyMixins(ColdObservable, [SubscriptionLoggable])`),
  // which silently breaks TestScheduler. The legacy webpack build was 'development'.
  const mode = env.NODE_ENV ?? 'development';

  const swcTs = {
    jsc: {
      target: 'es2015',
      parser: {syntax: 'typescript', tsx: false, decorators: true},
      transform: {legacyDecorator: true, decoratorMetadata: true},
    },
  };

  const includePaths = [path.resolve(__dirname, 'src'), ...SOURCE_LIBS];

  return {
    mode,
    cache: false,
    entry: {
      test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
      package: './src/package.ts',
    },
    experiments: {css: false},
    resolve: {
      extensions: ['.wasm', '.ts', '.mjs', '.js', '.json'],
      // Force every `rxjs` import to resolve to LibTests' own copy. Required so
      // TestScheduler can intercept AsyncScheduler.delegate (single rxjs instance).
      alias: {
        'rxjs': path.resolve(__dirname, 'node_modules/rxjs'),
      },
    },
    module: {
      rules: [
        {test: /\.ts$/, include: includePaths, loader: 'builtin:swc-loader', options: swcTs},
        {test: /\.css$/, use: ['style-loader', 'css-loader'], type: 'javascript/auto'},
        {test: /\.(mjs|js)$/, enforce: 'pre', use: ['source-map-loader']},
      ],
    },
    plugins: [
      new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
    ],
    devtool: 'source-map',
    // Note: `rxjs` is intentionally NOT externalized — TestScheduler needs to share
    // the same rxjs instance with code under test, so it must be bundled in.
    externals: {
      'datagrok-api/dg': 'DG',
      'DG': 'DG',
      'datagrok-api/grok': 'grok',
      'datagrok-api/ui': 'ui',
      'openchemlib/full.js': 'OCL',
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
    // splitChunks/runtimeChunk off: keep the fitting worker chunk self-contained.
    // A split vendor chunk (e.g. dayjs) would force the classic worker to
    // importScripts() it, which needs a publicPath the worker can't auto-detect
    // (no document) — the worker then fails to load. Webpack bundled dayjs into
    // the worker; this preserves that.
    optimization: {concatenateModules: false, splitChunks: false, runtimeChunk: false},
  };
};
