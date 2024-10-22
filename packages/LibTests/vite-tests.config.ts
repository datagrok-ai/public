// vite.config.js
import path, {resolve} from 'path';
import {defineConfig} from 'vite';
import cssInjectedByJsPlugin from 'vite-plugin-css-injected-by-js';

export default defineConfig({
  plugins: [
    cssInjectedByJsPlugin({
      jsAssetsFilterFunction: () => true,
    }),
  ],
  resolve: {
    alias: {
      'rxjs': path.resolve('node_modules/rxjs'),
      'rxjs/operator': path.resolve('node_modules/rxjs/operator'),
    },
    extensions: ['.wasm', '.mjs', '.ts', '.tsx', '.js', '.json'],
  },
  build: {
    emptyOutDir: false,
    sourcemap: true,
    minify: true,
    lib: {
      // Could also be a dictionary or array of multiple entry points
      entry: {
        'package-test': resolve(__dirname, 'src/package-test.ts'),
      },
      name: 'libtests_test',
      formats: ['iife'],
      fileName: (_, entryName) => `${entryName}.js`,
    },
    rollupOptions: {
      // make sure to externalize deps that shouldn't be bundled
      // into your library
      external: [
        'datagrok-api/dg',
        'DG',
        'datagrok-api/grok',
        'datagrok-api/ui',
        'openchemlib/full.js',
        'cash-dom',
        'dayjs',
        'wu',
        'exceljs',
        'html2canvas',
        'vue',
        'rxjs',
        'rxjs/operators',
        'rxjs/testing',
      ],
      output: {
        // Provide global variables to use in the UMD build
        // for externalized deps
        globals: {
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
          'rxjs': 'rxjs',
          'rxjs/operators': 'rxjs.operators',
          'rxjs/testing': 'rxjs.testing',
        },
        compact: true,
      },
    },
  },
});
