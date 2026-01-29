const path = require('path');
const { Compilation, sources } = require('webpack');

class AddAnnotationsPlugin {
  apply(compiler) {
    compiler.hooks.thisCompilation.tap('Replace', (compilation) => {
      compilation.hooks.processAssets.tap(
        {
          name: 'AddAnnotationsPlugin',
          stage: Compilation.PROCESS_ASSETS_STAGE_OPTIMIZE,
        },
        () => {
          const file = compilation.getAsset('preparser.js');
          const origSrc = file.source.source();
          const header = '//name: IvpPreparser\n//language: nodejs\n//input: string code\n//output: string result\n';
          compilation.updateAsset(
            'preparser.js',
            new sources.RawSource(header + origSrc.split("\n").slice(1).join("\n"))
          );
        }
      );
    });
  }
}

module.exports = {
  cache: {
    type: 'filesystem',
  },
  mode: 'production',
  entry: {
    main: './src/preparser.ts',
  },
  resolve: {
    extensions: ['.wasm', '.mjs', '.ts', '.json', '.js', '.tsx'],
  },
  plugins: [new AddAnnotationsPlugin()],
  module: {
    rules: [
      { test: /\.tsx?$/, loader: 'ts-loader', options: { allowTsInNodeModules: true } },
    ],
  },
  devtool: false,
  output: {
    filename: 'preparser.js',
    path: path.resolve(__dirname, 'scripts'),
    iife: false,
    environment: {
      templateLiteral: false,
    }
  },
  optimization: {
    minimize: false
  }
};
