const path = require('path');

module.exports = (env) => {
  return {
    mode: 'production',
    entry: {
      package: './src/index.ts',
    },
    resolve: {
      extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
    },
    module: {
      rules: [
	{test: /\.tsx?$/, loader: 'babel-loader', options: {
	  'plugins': ['@vue/babel-plugin-jsx']
	}},
	{test: /\.tsx?$/, loader: 'ts-loader'},
	{
          test: /\.(js|mjs|jsx|ts|tsx)$/,
          enforce: 'pre',
          use: ['source-map-loader'],
	}
      ],
    },
    devtool: 'source-map',
    externals: {
      'datagrok-api/dg': 'datagrok-api/dg',
      'datagrok-api/grok': 'datagrok-api/grok',
      'datagrok-api/ui': 'datagrok-api/ui',
      'openchemlib/full.js': 'openchemlib/full.js',
      'rxjs': 'rxjs',
      'rxjs/operators': 'rxjs/operators',
      'cash-dom': 'cash-dom',
      'dayjs': 'dayjs',
      'wu': 'wu',
      'exceljs': 'exceljs',
      'html2canvas': 'html2canvas',
    },
    experiments: {
      outputModule: true,
    },
    output: {
      path: path.resolve(__dirname, 'dist'),
      filename: 'index.mjs',
      library: {
	type: 'module'
      }
    },
    optimization: {
      concatenateModules: false,
      usedExports: true,
    }
  };
}
