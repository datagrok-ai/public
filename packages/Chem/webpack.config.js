const rdkitLibVersion = require('./src/rdkit_lib_version.js');
const path = require('path');
const {execSync} = require('child_process');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

function getDatagrokTools() {
  const pluginPath = 'datagrok-tools/plugins/func-gen-plugin';
  try {
    return require(pluginPath);
  } catch (e) {
    try {
      const globalPath = execSync('npm root -g').toString().trim();
      return require(path.join(globalPath, pluginPath));
    } catch (globalErr) {
      console.error('\n' + '='.repeat(60));
      console.error('ERROR: datagrok-tools not found!');
      console.error('To fix this, please install the tools globally by running:');
      console.error('\n    npm install -g datagrok-tools\n');
      console.error('='.repeat(60) + '\n');
      process.exit(1);
    }
  }
}

const FuncGeneratorPlugin = getDatagrokTools();

const mode = process.env.NODE_ENV ?? 'production';
if (mode !== 'production')
  console.warn(`Building Chem in '${mode}' mode.`);

const externalLibs = {
  'datagrok-api/dg': 'DG',
  'datagrok-api/grok': 'grok',
  'datagrok-api/ui': 'ui',
  'openchemlib/full': 'OCL',
  'rxjs': 'rxjs',
  'rxjs/operators': 'rxjs.operators',
  'cash-dom': '$',
  'wu': 'wu',
  'dayjs': 'dayjs',
  'NGL': 'NGL',
};

module.exports = (_env, _options) => ({
  cache: {type: 'filesystem'},
  stats: {
    children: true,
  },
  mode: mode,
  entry: {
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
    package: [`./src/${rdkitLibVersion}.wasm`, './src/bbwasm_bg.wasm', './src/package.ts'],
  },
  devtool: 'source-map',
  devServer: {
    contentBase: './dist',
  },
  externals: [
    function({context, request}, callback) {
      if (context.includes('worker') && request.includes('openchemlib/full')) {
        //console.log({context, request});
        return callback();
      }

      if (externalLibs[request])
        return callback(null, externalLibs[request], 'var');
      return callback();
    },
  ],
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
  resolve: {
    fallback: {
      url: false,
      perf_hooks: false,
    },
    extensions: ['.wasm', '.ts', '.mjs', '.js', '.json', '.tsx'],
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
  module: {
    rules: [
      // prevent a webpack name change
      {
        test: /\.(wasm)$/i,
        type: 'javascript/auto',
        loader: 'file-loader',
        options: {
          publicPath: 'dist/',
          name: '[name].[ext]',
        },
      },
      {
        test: /\.ts$/,
        loader: 'ts-loader',
        exclude: /node_modules/,
      },
      {
        test: /\.js$/,
        enforce: 'pre',
        use: ['source-map-loader'],
        exclude: /node_modules/,
      },
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader'],
      },
    ],
  },
});

/**
 * hacky way so that the externals are correctly detected by script in grok check
 externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'wu': 'wu',
    'dayjs': 'dayjs',
    'NGL': 'NGL',
}
 */
