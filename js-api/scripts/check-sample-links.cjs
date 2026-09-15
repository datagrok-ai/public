// Checks every ApiSamples link in the JSDoc against packages/ApiSamples/scripts.
// A link must name an existing sample; exits 1 otherwise. Run: npm run check-links
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const samplesRoot = path.resolve(root, '../packages/ApiSamples/scripts');
if (!fs.existsSync(samplesRoot)) {
  console.log(`check-links: ${samplesRoot} not checked out, skipping`);
  process.exit(0);
}

const files = [path.join(root, 'ui.ts')];
(function walk(d) {
  for (const e of fs.readdirSync(d, {withFileTypes: true})) {
    const p = path.join(d, e.name);
    if (e.isDirectory()) { if (!/node_modules|[\\/]datagrok$|[\\/]node$|interfaces/.test(p)) walk(p); }
    else if (/\.ts$/.test(e.name) && !/\.d\.ts$|\.g\.ts$/.test(e.name)) files.push(p);
  }
})(path.join(root, 'src'));

const samples = new Set();
(function walk(d) {
  for (const e of fs.readdirSync(d, {withFileTypes: true})) {
    const p = path.join(d, e.name);
    if (e.isDirectory()) walk(p);
    else samples.add(path.relative(samplesRoot, p).replace(/\\/g, '/').replace(/\.(js|ts|py|r|jl|md)$/, ''));
  }
})(samplesRoot);

const re = /https?:\/\/(?:public|dev)\.datagrok\.ai\/(?:js\/samples|script\/samples\/javascript)\/([a-zA-Z0-9/_.\- ]+?)(?:[)\s}|]|$)/g;
let total = 0;
const dead = [];
for (const f of files) {
  fs.readFileSync(f, 'utf8').split('\n').forEach((line, i) => {
    let m;
    while ((m = re.exec(line)) !== null) {
      total++;
      const target = m[1].trim();
      if (!samples.has(target))
        dead.push(`${path.relative(root, f).replace(/\\/g, '/')}:${i + 1}  ${target}`);
    }
  });
}
console.log(`check-links: ${total} sample links, ${dead.length} broken`);
for (const d of dead) console.log('  ' + d);
process.exit(dead.length ? 1 : 0);
