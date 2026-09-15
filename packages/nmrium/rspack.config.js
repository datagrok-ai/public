const path = require('path');
const {bundler} = require('@datagrok/build-config');

module.exports = bundler({
  jsx: 'react',
  resolve: {
    alias: {
      'openchemlib/full$': 'openchemlib/full.js',
      './lib/components/StructureEditor': './lib/components/StructureEditor.js',
      './lib/components/SmilesSvgRenderer': './lib/components/SmilesSvgRenderer.js',
      './lib/components/MolfileSvgRenderer': './lib/components/MolfileSvgRenderer.js',
      './lib/components/IdcodeSvgRenderer': './lib/components/IdcodeSvgRenderer.js',
      'cheminfo-font/lib-react-cjs/lib-react-tsx/nmr/Peaks': path.resolve(__dirname, 'node_modules/cheminfo-font/lib-react-cjs/lib-react-tsx/nmr/Peaks.js'),
    },
  },
});
