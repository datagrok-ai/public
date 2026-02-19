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
    test: {
      filename: 'package-test.js',
      library: {type: 'var', name: `${packageName}_test`},
      import: './src/package-test.ts',
    },
  },
  devtool: 'inline-source-map',
  devServer: {
    contentBase: './dist',
  },
  target: 'web',
  module: {
    rules: [
      {
        test: /\.ts?$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
    ],
  },
  resolve: {
    extensions: ['.mjs', '.ts', '.json', '.js', '.tsx'],
  },
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full': 'OCL',
  },
  output: {
    filename: '[name].js',
    library: 'drugbank',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins: [
    new FuncGeneratorPlugin({outputPath: './src/package.g.ts'}),
  ],
};
