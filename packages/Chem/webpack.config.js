const path = require('path');
const rdkitLibVersion = require('./src/rdkit_lib_version.js');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

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
    package: [`./src/${rdkitLibVersion}.wasm`, './src/package.ts'],
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
        test: /\.worker\.ts$/,
        loader: 'worker-loader',
        options: {
          inline: 'fallback', // this creates a separate file
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
