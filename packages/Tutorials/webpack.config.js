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
  
module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'development',
  entry: {
    package: './src/package.ts',
    test: {filename: 'package-test.js', library: {type: 'var', name:`${packageName}_test`}, import: './src/package-test.ts'},
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
  },
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader' },
      { test: /\.css$/i, use: ['style-loader', 'css-loader'] },
    ],
  },
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
  },
  output: {
    filename: '[name].js',
    library: 'tutorials',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins:[
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'})
  ]
};
