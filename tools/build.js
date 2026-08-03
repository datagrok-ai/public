#!/usr/bin/env node
// Transpiles bin/**/*.ts -> bin/**/*.js in place using @babel/core.
// Replaces `@babel/cli` (which pulls deprecated glob@7/inflight); the babel
// presets/plugins are read from the "babel" key in package.json automatically.
const fs = require('fs');
const path = require('path');
const babel = require('@babel/core');

const binDir = path.join(__dirname, 'bin');
const sourceMaps = process.argv.includes('--source-maps');

function collectTsFiles(dir, acc) {
  for (const entry of fs.readdirSync(dir, {withFileTypes: true})) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      // Tests run from .ts via vitest and are never required at runtime.
      if (entry.name === '__tests__')
        continue;
      collectTsFiles(full, acc);
    } else if (entry.isFile() && entry.name.endsWith('.ts'))
      acc.push(full);
  }
  return acc;
}

let count = 0;
for (const file of collectTsFiles(binDir, [])) {
  const result = babel.transformFileSync(file, {sourceMaps: sourceMaps});
  if (!result || result.code == null)
    continue;
  const outFile = file.replace(/\.ts$/, '.js');
  let code = result.code;
  if (sourceMaps && result.map) {
    const mapFile = outFile + '.map';
    fs.writeFileSync(mapFile, JSON.stringify(result.map));
    code += `\n//# sourceMappingURL=${path.basename(mapFile)}\n`;
  }
  fs.writeFileSync(outFile, code);
  count++;
}
console.log(`Successfully compiled ${count} files with Babel.`);
