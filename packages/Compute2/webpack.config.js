const path = require('path');
const packageName = path.parse(require('./package.json').name).name.toLowerCase().replace(/-/g, '');
const webpack = require('webpack');

module.exports = (env) => {
  const ENABLE_VUE_DEV_TOOLS = env.enable_vue_dev_tools;
  const mode = ENABLE_VUE_DEV_TOOLS ? 'development' : (env.NODE_ENV ?? 'production');

  if (ENABLE_VUE_DEV_TOOLS)
    console.warn('Building DEV ONLY Compute2 build with vue devtools support');

  return {
    mode,
    entry: {
      test: {filename: 'package-test.js', library: {type: 'var', name: `${packageName}_test`}, import: './src/package-test.ts'},
      package: './src/package.ts',
    },
    resolve: {
      extensions: ['.wasm', '.ts', '.tsx', '.mjs', '.js', '.jsx',  '.json'],
      alias: {
	...(ENABLE_VUE_DEV_TOOLS ? {'vue': path.resolve('node_modules/vue/dist/vue.esm-bundler.js')} : {}),
      }
    },
    module: {
      rules: [
	{test: /\.tsx?$/, loader: 'babel-loader', options: {
          'plugins': [['@vue/babel-plugin-jsx', { isCustomElement: tag => tag.startsWith('dg-') || tag === 'dock-spawn-ts' }]],
	}},
	{test: /\.tsx?$/, loader: 'ts-loader', options: {allowTsInNodeModules: true}},
	{
          test: /\.css$/,
          use: ['style-loader', 'css-loader', 'postcss-loader'],
	},
	{
          test: /\.(ts|tsx|mjs|js|jsx)$/,
          enforce: 'pre',
          use: ['source-map-loader'],
	},
      ],
    },
    plugins: [new webpack.DefinePlugin({
      ENABLE_VUE_DEV_TOOLS,
      __VUE_PROD_HYDRATION_MISMATCH_DETAILS__: 'false',
      __VUE_PROD_DEVTOOLS__: 'false',
      __VUE_OPTIONS_API__: 'true',
    })],
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
      ...(ENABLE_VUE_DEV_TOOLS ? {} : {'vue': 'Vue'}),
    },
    output: {
      filename: '[name].js',
      library: packageName,
      libraryTarget: 'var',
      path: path.resolve(__dirname, 'dist'),
      clean: true,
    },
    optimization: {
      usedExports: true,
      concatenateModules: false,
    }
  }
};
