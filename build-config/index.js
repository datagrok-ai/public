// Shared rspack configuration for Datagrok plugins.
//
//   // rspack.config.js (only needed when a package deviates from the defaults)
//   module.exports = require('@datagrok/build-config').bundler({externals: {'ngl': 'NGL'}});
//
// A package without rspack.config.js gets bundler({}). Every loader is resolved from this
// package so plugins carry no bundler dependencies of their own.
const fs = require('fs');
const path = require('path');
const {rspack} = require('@rspack/core');

// Provided by the platform at runtime; never bundled.
const PLATFORM_EXTERNALS = {
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
};

const ASSET_TEST = /\.(png|jpe?g|gif|svg|ico|sdf|mol|woff2?|ttf|eot|otf|csv|txt|md)$/;

function swcOptions({jsx, decorators = true}) {
  return {
    jsc: {
      target: 'es2015',
      parser: {syntax: 'typescript', decorators, tsx: jsx === 'react'},
      transform: {
        legacyDecorator: decorators,
        decoratorMetadata: decorators,
        useDefineForClassFields: false,
        react: jsx === 'react' ? {runtime: 'automatic'} : undefined,
      },
    },
  };
}

function mergeConfig(base, o) {
  const out = {...base};
  for (const [k, v] of Object.entries(o)) {
    if (['dir', 'jsx', 'wasm', 'decorators', 'assets', 'rules', 'mode', 'css'].includes(k)) continue;
    if (k === 'externals') {
      out.externals = {...base.externals, ...v};
      for (const [name, value] of Object.entries(v)) if (value === false) delete out.externals[name];
    }
    else if (k === 'resolve') out.resolve = {...base.resolve, ...v, alias: {...base.resolve.alias, ...(v.alias || {})}, fallback: {...base.resolve.fallback, ...(v.fallback || {})}};
    else if (k === 'plugins') out.plugins = [...base.plugins, ...v];
    else if (k === 'experiments') out.experiments = {...base.experiments, ...v};
    else if (k === 'output') out.output = {...base.output, ...v};
    else if (k === 'optimization') out.optimization = {...base.optimization, ...v};
    else if (k === 'module') out.module = {...base.module, ...v, rules: [...base.module.rules, ...((v && v.rules) || [])]};
    else out[k] = v;
  }
  return out;
}

/**
 * Build an rspack config for the package in `o.dir` (default: cwd).
 * Options beyond rspack's own: `jsx: 'react'`, `wasm: 'async' | 'sync' | 'asset'`,
 * `decorators: false`, `assets: RegExp` (extra asset/resource test), `rules: []` (prepended).
 */
function bundler(o = {}) {
  const dir = o.dir || process.cwd();
  const pkg = JSON.parse(fs.readFileSync(path.join(dir, 'package.json'), 'utf8'));
  const name = path.parse(pkg.name).name.toLowerCase().replace(/-/g, '');
  const mode = o.mode || (process.env.NODE_ENV === 'development' ? 'development' : 'production');
  const hasTests = fs.existsSync(path.join(dir, 'src', 'package-test.ts'));
  const {FuncGeneratorPlugin} = loadFuncGen();

  const ext = fs.existsSync(path.join(dir, 'src', 'package.ts')) ? 'ts' : 'js';
  const entry = {package: `./src/package.${ext}`};
  const testEntry = ['ts', 'js'].map((e) => `./src/package-test.${e}`).find((e) => fs.existsSync(path.join(dir, e)));
  if (testEntry)
    entry.test = {filename: 'package-test.js', library: {type: 'var', name: `${name}_test`}, import: testEntry};

  const rules = [
    ...(o.rules || []),
    {test: /\.tsx?$/, exclude: /node_modules/, loader: 'builtin:swc-loader', options: swcOptions(o)},
    // Library dist is ES2020; downlevel it for the platform's oldest browser like the .ts sources.
    {test: /[\\/]libraries[\\/][^\\/]+[\\/]dist[\\/].*\.js$/, loader: 'builtin:swc-loader',
      options: {env: {targets: 'chrome 50'}, jsc: {parser: {syntax: 'ecmascript'}}}},
    ...(o.css === false ? [] : [{test: /\.css$/, use: [require.resolve('style-loader'), require.resolve('css-loader')], type: 'javascript/auto'}]),
    {test: o.assets ? new RegExp(`${ASSET_TEST.source}|${o.assets.source}`) : ASSET_TEST, type: 'asset/resource'},
  ];
  if (o.wasm === 'asset' || !o.wasm)
    rules.push({test: /\.wasm$/, type: 'asset/resource'});

  const base = {
    context: dir,
    mode,
    entry,
    devtool: 'source-map',
    resolve: {
      extensions: ['.ts', '.tsx', '.mjs', '.js', '.jsx', '.json', '.wasm'],
      // `import './x.worker.ts'` inside a compiled library resolves to the emitted x.worker.js
      extensionAlias: {'.ts': ['.ts', '.js'], '.js': ['.js', '.ts']},
      alias: {},
      fallback: {url: false, fs: false, path: false, crypto: false},
    },
    module: {rules},
    externals: {...PLATFORM_EXTERNALS},
    plugins: [
      new rspack.DefinePlugin({'process.env.NODE_ENV': JSON.stringify(mode), 'process.env': '{}'}),
      ...(FuncGeneratorPlugin ? [new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})] : []),
    ],
    output: {
      filename: '[name].js',
      library: {name, type: 'var'},
      path: path.resolve(dir, 'dist'),
      clean: true,
    },
    optimization: {minimize: mode === 'production'},
    experiments: {
      css: false,
      topLevelAwait: true,
      asyncWebAssembly: o.wasm === 'async',
      syncWebAssembly: o.wasm === 'sync',
    },
    stats: 'errors-warnings',
  };
  return mergeConfig(base, o);
}

// The plugin generates src/package.g.ts and src/package-api.ts on every bundle; a build without it
// would silently ship stale function metadata, so failing to load it is fatal.
function loadFuncGen() {
  try {
    return {FuncGeneratorPlugin: require('datagrok-tools/plugins/func-gen-plugin')};
  }
  catch (e) {
    throw new Error(`[build-config] cannot load datagrok-tools/plugins/func-gen-plugin: ${e.message}`);
  }
}

// Absolute loader paths for hand-written configs: plugins carry no loader dependencies.
const loaders = {
  style: require.resolve('style-loader'),
  css: require.resolve('css-loader'),
  null: require.resolve('null-loader'),
};

module.exports = {bundler, PLATFORM_EXTERNALS, swcOptions, rspack, loaders};
