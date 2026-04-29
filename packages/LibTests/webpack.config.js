const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');

module.exports = {
  cache: {
    type: 'filesystem',
    buildDependencies: {
      config: [__filename],
    },
  },
  mode: 'development',
  entry: {
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
    package: './src/package.ts',
  },
  resolve: {
    extensions: ['.wasm', '.ts', '.tsx', '.mjs', '.js', '.jsx',  '.json'],
    // Needed for rxjs test scheduler to work properly
    alias: {
      'rxjs': path.resolve('node_modules/rxjs'),
      'rxjs/operator': path.resolve('node_modules/rxjs/operator'),
    },
  },
  module: {
    rules: [
      {test: /\.tsx?$/, loader: 'ts-loader', options: {allowTsInNodeModules: true}},
      {
        test: /\.css$/,
        use: ['style-loader', 'css-loader'],
      },
      {
        test: /\.(ts|tsx|mjs|js|jsx)$/,
        enforce: 'pre',
        use: ['source-map-loader'],
      },
    ],
  },
  devtool: 'source-map',
  externals: [
    // Workers don't get the platform's globals (no `window.dayjs`), so the
    // fitting worker bundle (and the worker-safe `webworkers/` modules it
    // pulls in) imports `dayjs` and webpack must bundle it. The main
    // bundle still externalizes dayjs to the platform global.
    ({context, request}, cb) => {
      if (request === 'dayjs' && context &&
          /[\\/](?:fitting[\\/]worker|webworkers)(?:[\\/]|$)/.test(context))
        return cb();
      const platformExternals = {
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
        'vue': 'Vue',
      };
      if (request in platformExternals) return cb(null, platformExternals[request]);
      cb();
    },
  ],
  output: {
    filename: '[name].js',
    library: packageName,
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
    clean: true,
  },
};
