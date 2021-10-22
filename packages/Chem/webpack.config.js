const path = require('path');

module.exports = {
  mode: 'development',
  entry: {
    package: ['./src/RDKit_minimal_2021.03_17.wasm', './src/package.js'/*, './src/rdkit_worker.js'*/]
  },
  // devtool: 'inline-source-map',
  devServer: {
    contentBase: './dist'
  },
  externals: {
    'datagrok-api/dg': 'DG',
    'datagrok-api/grok': 'grok',
    'datagrok-api/ui': 'ui',
    'openchemlib/full.js': 'OCL',
    'rxjs': 'rxjs',
    'rxjs/operators': 'rxjs.operators',
    'cash-dom': '$',
    'dayjs': 'dayjs',
  },
  output: {
    filename: '[name].js',
    library: 'chem',
    libraryTarget: 'var',
    path: path.resolve(__dirname, 'dist'),
  },
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true
  },
  /*
  resolve: {
    extensions: ['.wasm', '.mjs', '.js', '.json', '.ts', '.tsx'],
  },
  */
  module: {
    rules: [
      // WASM files should not be processed but just
      // be emitted and we want to have their public URL.
      // file-loader -> options -> name =  '[name].[ext]'  - prevents webpack name change
      {
        test: /\.(wasm)$/i,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/",
          name: '[name].[ext]'
        }
      }/*,
      {
        test: /\_(worker.js)$/i,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/",
          name: '[name].[ext]'
        }
      }*/,
      {
        test: /_worker\.js$/i,
        loader: "worker-loader",
        options: {
          inline:"no-fallback"
        },
      },
      {
        test: /\.tsx?$/,
        loader: 'ts-loader'
      }
    ],
  }/*
  resolve: {
    fallback: {
      path: "false",
      fs: "false"
    },
  }*/
};
