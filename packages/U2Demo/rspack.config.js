const {bundler} = require('@datagrok/build-config');

const config = bundler({externals: {'datagrok-api/u2core': 'DG.U2'}});

// `?raw` inlines a source file as text so the demo can show the source of the build that is
// running (src/source-panel.ts); the TypeScript rule has to skip those requests, or it compiles
// them instead and the import resolves to a module with no default export.
for (const rule of config.module.rules)
  if (String(rule.test) === String(/\.tsx?$/))
    rule.resourceQuery = {not: [/raw/]};
config.module.rules.unshift({resourceQuery: /raw/, type: 'asset/source'});

module.exports = config;
