const path = require('path');
const fs = require('fs');
const rspack = require('@rspack/core');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const FuncGeneratorPlugin = require('datagrok-tools/plugins/func-gen-plugin');

// Resolve linked-or-npm-installed lib paths to their real on-disk location so
// module.rules[].include matches via the resolved path (not the node_modules symlink).
// Skips libs that aren't installed (e.g., transitive lib was uninstalled).
const realLib = (pkg) => {
  try { return fs.realpathSync(path.resolve(__dirname, 'node_modules', pkg)); }
  catch { return null; }
};
const SOURCE_LIBS = [
  '@datagrok-libraries/compute-utils',
  '@datagrok-libraries/webcomponents',
  '@datagrok-libraries/webcomponents-vue',
  '@datagrok-libraries/utils',
  '@datagrok-libraries/arrow',
  '@datagrok-libraries/dock-spawn-dg',
  '@datagrok-libraries/test',
  'diff-grok', // ships .ts source; transitive dep of compute-utils
].map(realLib).filter(Boolean);

module.exports = (env = {}) => {
  const ENABLE_VUE_DEV_TOOLS = env.enable_vue_dev_tools;
  const mode = ENABLE_VUE_DEV_TOOLS ? 'development' : (env.NODE_ENV ?? 'production');

  if (ENABLE_VUE_DEV_TOOLS)
    console.warn('Building DEV ONLY Compute2 build with vue devtools support');

  const swcCommon = {
    jsc: {
      target: 'es2015',
      parser: {syntax: 'typescript', decorators: true},
      transform: {legacyDecorator: true, decoratorMetadata: true},
    },
  };

  const swcTs = {
    ...swcCommon,
    jsc: {...swcCommon.jsc, parser: {...swcCommon.jsc.parser, tsx: false}},
  };

  const swcTsx = {
    ...swcCommon,
    jsc: {
      ...swcCommon.jsc,
      parser: {...swcCommon.jsc.parser, tsx: true},
      experimental: {
        plugins: [['swc-plugin-vue-jsx', {customElementPatterns: ['^dg-', '^dock-spawn-ts$']}]],
      },
    },
  };

  const swcJsx = {
    jsc: {
      target: 'es2015',
      parser: {syntax: 'ecmascript', jsx: true},
      experimental: {
        plugins: [['swc-plugin-vue-jsx', {customElementPatterns: ['^dg-', '^dock-spawn-ts$']}]],
      },
    },
  };

  const includePaths = [path.resolve(__dirname, 'src'), ...SOURCE_LIBS];

  const config = {
    mode,
    cache: false,
    entry: {
      test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
      package: './src/package.ts',
    },
    experiments: {css: false},
    resolve: {
      // .ts first → deep imports into linked libs resolve to source; npm deps still fall through to .js because their .ts doesn't exist
      extensions: ['.wasm', '.ts', '.tsx', '.mjs', '.js', '.jsx', '.json'],
      // 'source' field comes first → top-level package imports go to ./index.ts
      mainFields: ['source', 'browser', 'module', 'main'],
      alias: {
        ...(ENABLE_VUE_DEV_TOOLS ? {'vue': path.resolve('node_modules/vue/dist/vue.esm-bundler.js')} : {}),
      },
    },
    module: {
      rules: [
        {test: /\.tsx$/, include: includePaths, loader: 'builtin:swc-loader', options: swcTsx},
        {test: /\.ts$/, include: includePaths, loader: 'builtin:swc-loader', options: swcTs},
        {test: /\.jsx$/, loader: 'builtin:swc-loader', options: swcJsx},
        {test: /\.css$/, use: ['style-loader', 'css-loader', 'postcss-loader'], type: 'javascript/auto'},
        {test: /\.(mjs|js)$/, enforce: 'pre', use: ['source-map-loader']},
      ],
    },
    plugins: [
      new rspack.DefinePlugin({
        ENABLE_VUE_DEV_TOOLS,
        __VUE_PROD_HYDRATION_MISMATCH_DETAILS__: 'false',
        __VUE_PROD_DEVTOOLS__: 'false',
        __VUE_OPTIONS_API__: 'true',
      }),
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
      'vue': 'Vue',
    },
    output: {
      filename: '[name].js',
      library: packageName,
      libraryTarget: 'var',
      path: path.resolve(__dirname, 'dist'),
      clean: true,
    },
    optimization: {concatenateModules: false},
  };

  if (ENABLE_VUE_DEV_TOOLS)
    delete config.externals.vue;

  return config;
};
