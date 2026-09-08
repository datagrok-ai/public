// Inventory of the datagrok-api public surface: how much of it is documented, tested (ApiTests),
// sampled (ApiSamples) and mentioned in help/develop. Name-based, so the coverage figures are upper bounds.
//
//   node scripts/inventory.cjs                    print the summary
//   node scripts/inventory.cjs --check            also fail when a ratchet metric got worse than inventory-baseline.json
//   node scripts/inventory.cjs --update-baseline  rewrite the baseline from the current numbers
//   node scripts/inventory.cjs --out <dir>        also write the tables (members.json, file-stats.md, class-stats.md,
//                                                 dark-members.txt, samples-by-member.json) into <dir>
//
// `npm run inventory` is the --check form and runs in the JS API workflow; `grok-core audit js-api` wraps it.
const ts = require('typescript');
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const args = process.argv.slice(2);
const check = args.includes('--check');
const updateBaseline = args.includes('--update-baseline');
const outIx = args.indexOf('--out');
const outDir = outIx >= 0 ? path.resolve(args[outIx + 1]) : null;
const baselineFile = path.join(__dirname, 'inventory-baseline.json');

const files = [];
(function walk(d) {
  for (const e of fs.readdirSync(d, {withFileTypes: true})) {
    const p = path.join(d, e.name);
    if (e.isDirectory()) { if (!/node_modules|[\\/]datagrok$|[\\/]node$/.test(p)) walk(p); }
    else if (/\.ts$/.test(e.name) && !/\.d\.ts$/.test(e.name)) files.push(p);
  }
})(path.join(root, 'src'));
for (const f of ['ui.ts', 'grok.ts', 'dg.ts']) files.push(path.join(root, f));

function readAll(d) {
  let s = '';
  (function walk(dd) {
    for (const e of fs.readdirSync(dd, {withFileTypes: true})) {
      const p = path.join(dd, e.name);
      if (e.isDirectory()) { if (!/node_modules|dist/.test(e.name)) walk(p); }
      else if (/\.(ts|js|md|mdx)$/.test(e.name)) s += fs.readFileSync(p, 'utf8') + '\n';
    }
  })(d);
  return s;
}
const testsText = readAll(path.resolve(root, '../packages/ApiTests/src'));
const samplesText = readAll(path.resolve(root, '../packages/ApiSamples/scripts'));
const helpText = readAll(path.resolve(root, '../help/develop'));

const members = [];
const fileStats = {};

function hasJsDoc(node) { return ts.getJSDocCommentsAndTags(node).length > 0; }
function docText(node) { return ts.getJSDocCommentsAndTags(node).map((t) => t.getFullText()).join('\n'); }
function isExported(node) { return !!(ts.getCombinedModifierFlags(node) & ts.ModifierFlags.Export); }
function isPublic(node) {
  const f = ts.getCombinedModifierFlags(node);
  return !(f & ts.ModifierFlags.Private) && !(f & ts.ModifierFlags.Protected);
}
function nameOf(n) { return n.name ? n.name.getText() : '(anon)'; }

for (const file of files) {
  const src = fs.readFileSync(file, 'utf8');
  const sf = ts.createSourceFile(file, src, ts.ScriptTarget.ES2020, true);
  const rel = path.relative(root, file).replace(/\\/g, '/');
  const st = fileStats[rel] = {
    loc: src.split('\n').length, exports: 0, publicMembers: 0, documented: 0,
    any: (src.match(/\bany\b/g) || []).length,
    tsIgnore: (src.match(/@ts-ignore/g) || []).length,
    todo: (src.match(/TODO|FIXME|HACK/g) || []).length,
    deprecated: (src.match(/@deprecated/gi) || []).length,
    eslintDisable: (src.match(/eslint-disable/g) || []).length,
  };

  function addMember(owner, kind, node, name) {
    const doc = hasJsDoc(node);
    const sig = src.slice(node.getStart(), Math.min(node.getEnd(), node.getStart() + 160)).split('{')[0].replace(/\s+/g, ' ').trim();
    const dt = docText(node);
    members.push({file: rel, owner, kind, name, hasDoc: doc, docLen: dt.length, deprecated: /@deprecated/i.test(dt), anyInSig: /\bany\b/.test(sig), sig});
    st.publicMembers++;
    if (doc) st.documented++;
  }

  ts.forEachChild(sf, function visit(node) {
    if (ts.isClassDeclaration(node) || ts.isInterfaceDeclaration(node)) {
      if (!isExported(node)) return;
      st.exports++;
      const owner = nameOf(node);
      const isIface = ts.isInterfaceDeclaration(node);
      addMember('', isIface ? 'interface' : 'class', node, owner);
      for (const m of node.members) {
        if (!isPublic(m)) continue;
        if (ts.isConstructorDeclaration(m)) continue;
        const nm = m.name ? m.name.getText() : '(index)';
        if (nm.startsWith('_') || nm === 'dart') continue;
        const kind = ts.isMethodDeclaration(m) || ts.isMethodSignature(m) ? 'method' : ts.isGetAccessor(m) ? 'getter' : ts.isSetAccessor(m) ? 'setter' : 'prop';
        if (kind === 'setter') continue;
        addMember(owner, kind, m, nm);
      }
    }
    else if (ts.isFunctionDeclaration(node)) {
      if (!isExported(node)) return;
      st.exports++;
      addMember('', 'function', node, nameOf(node));
    }
    else if (ts.isEnumDeclaration(node) || ts.isTypeAliasDeclaration(node)) {
      if (!isExported(node)) return;
      st.exports++;
      addMember('', ts.isEnumDeclaration(node) ? 'enum' : 'type', node, nameOf(node));
    }
    else if (ts.isVariableStatement(node)) {
      if (!isExported(node)) return;
      for (const d of node.declarationList.declarations) {
        st.exports++;
        addMember('', 'const', node, d.name.getText());
      }
    }
    else if (ts.isModuleDeclaration(node)) {
      if (node.body) ts.forEachChild(node.body, visit);
    }
  });
}

function esc(s) { return s.replace(/[$.]/g, '\\$&'); }
function covered(text, m) {
  if (m.owner === '')
    return new RegExp('\\b' + esc(m.name) + '\\b').test(text);
  return new RegExp('\\.' + esc(m.name) + '\\b').test(text);
}
for (const m of members) {
  m.inTests = covered(testsText, m);
  m.inSamples = covered(samplesText, m);
  m.inHelp = covered(helpText, m);
  m.generated = /\.g\.ts$|src\/interfaces\//.test(m.file);
}

// Samples per member: an `//api: DG.Class.member, ui.fn, grok.ns.fn` header line names the members a
// sample demonstrates; without one, qualified usages (`DG.Viewer.fromType`, `ui.input.string`, `grok.data.demo`)
// are taken from the source. Bare `.member(` calls are ambiguous and not counted.
const samplesDir = path.resolve(root, '../packages/ApiSamples/scripts');
const samplesByMember = {};
let headed = 0;
(function walk(dd) {
  for (const e of fs.readdirSync(dd, {withFileTypes: true})) {
    const p = path.join(dd, e.name);
    if (e.isDirectory()) walk(p);
    else if (/\.js$/.test(e.name)) {
      const text = fs.readFileSync(p, 'utf8');
      const rel = path.relative(samplesDir, p).replace(/\\/g, '/').replace(/\.js$/, '');
      const header = text.match(/^\/\/\s*api:\s*(.+)$/m);
      if (header) headed++;
      const names = header
        ? header[1].split(',').map((s) => s.trim()).filter(Boolean)
        : [...new Set([...text.matchAll(/\b((?:DG|ui|grok)(?:\.[A-Za-z_$][\w$]*){1,3})/g)].map((x) => x[1]))];
      for (const n of names) (samplesByMember[n] = samplesByMember[n] || []).push(rel);
    }
  }
})(samplesDir);

const pct = (a, b) => (100 * a / b).toFixed(0) + '%';
const isMember = (m) => m.kind === 'method' || m.kind === 'getter' || m.kind === 'prop';
const dark = (m) => m.owner && !m.hasDoc && !m.inTests && !m.inSamples;
const hand = members.filter((m) => !m.generated);
const handMembers = hand.filter(isMember);
const handClasses = hand.filter((m) => m.kind === 'class');
const allMembers = members.filter(isMember);
const allClasses = members.filter((m) => m.kind === 'class');

// Ratchet metrics: counts that must not grow (a new member must arrive documented, typed and covered).
const metrics = {
  undocumentedMembers: handMembers.filter((m) => !m.hasDoc).length,
  undocumentedClasses: handClasses.filter((m) => !m.hasDoc).length,
  darkMembers: hand.filter(dark).length,
  anyInSignature: handMembers.filter((m) => m.anyInSig).length,
  undocumentedGenerated: members.filter((m) => m.generated && !m.hasDoc).length,
};

console.log('--- hand-written only ---');
console.log('entries', hand.length, 'documented', pct(hand.filter((m) => m.hasDoc).length, hand.length));
console.log('classes', handClasses.length, 'members', handMembers.length, 'documented', pct(handMembers.filter((m) => m.hasDoc).length, handMembers.length));
console.log('members in tests', pct(handMembers.filter((m) => m.inTests).length, handMembers.length), 'in samples', pct(handMembers.filter((m) => m.inSamples).length, handMembers.length), 'in help', pct(handMembers.filter((m) => m.inHelp).length, handMembers.length));
console.log('classes in tests', pct(handClasses.filter((m) => m.inTests).length, handClasses.length), 'in samples', pct(handClasses.filter((m) => m.inSamples).length, handClasses.length));
console.log('--- all incl generated ---');
console.log(`files ${files.length}, surface entries ${members.length}, documented ${pct(members.filter((m) => m.hasDoc).length, members.length)}`);
console.log(`classes ${allClasses.length}, members ${allMembers.length}, members documented ${pct(allMembers.filter((m) => m.hasDoc).length, allMembers.length)}`);
console.log(`members in tests ${pct(allMembers.filter((m) => m.inTests).length, allMembers.length)}, in samples ${pct(allMembers.filter((m) => m.inSamples).length, allMembers.length)}, in help ${pct(allMembers.filter((m) => m.inHelp).length, allMembers.length)}`);
console.log(`deprecated ${members.filter((m) => m.deprecated).length}; samples with an //api: header ${headed}; qualified names with a sample ${Object.keys(samplesByMember).length}`);
console.log('--- ratchet ---');
console.log(Object.entries(metrics).map(([k, v]) => `${k} ${v}`).join(', '));

if (outDir) {
  fs.mkdirSync(outDir, {recursive: true});
  fs.writeFileSync(path.join(outDir, 'members.json'), JSON.stringify(members));
  fs.writeFileSync(path.join(outDir, 'file-stats.json'), JSON.stringify(fileStats, null, 1));
  fs.writeFileSync(path.join(outDir, 'samples-by-member.json'), JSON.stringify(samplesByMember, null, 1));

  let md = '| file | LOC | exports | public members | doc % | any | ts-ignore | TODO | deprecated |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n';
  for (const [f, s] of Object.entries(fileStats).sort((a, b) => b[1].publicMembers - a[1].publicMembers))
    md += `| ${f} | ${s.loc} | ${s.exports} | ${s.publicMembers} | ${s.publicMembers ? (100 * s.documented / s.publicMembers).toFixed(0) : '-'} | ${s.any} | ${s.tsIgnore} | ${s.todo} | ${s.deprecated} |\n`;
  fs.writeFileSync(path.join(outDir, 'file-stats.md'), md);

  const byOwner = {};
  for (const m of members) {
    if (m.kind === 'class' || m.kind === 'interface') {
      byOwner[m.name] = Object.assign(byOwner[m.name] || {n: 0, doc: 0, tests: 0, samples: 0, help: 0}, {file: m.file, kind: m.kind, classDoc: m.hasDoc, classInTests: m.inTests, classInSamples: m.inSamples});
      continue;
    }
    if (!m.owner) continue;
    const o = byOwner[m.owner] = byOwner[m.owner] || {file: m.file, kind: '?', n: 0, doc: 0, tests: 0, samples: 0, help: 0};
    o.n++; if (m.hasDoc) o.doc++; if (m.inTests) o.tests++; if (m.inSamples) o.samples++; if (m.inHelp) o.help++;
  }
  let cmd = '| class | file | members | doc % | in tests % | in samples % | class in tests | class in samples |\n|---|---|---:|---:|---:|---:|:-:|:-:|\n';
  for (const [c, o] of Object.entries(byOwner).sort((a, b) => b[1].n - a[1].n))
    cmd += `| ${c} | ${o.file} | ${o.n} | ${o.n ? (100 * o.doc / o.n).toFixed(0) : '-'} | ${o.n ? (100 * o.tests / o.n).toFixed(0) : '-'} | ${o.n ? (100 * o.samples / o.n).toFixed(0) : '-'} | ${o.classInTests ? 'y' : ''} | ${o.classInSamples ? 'y' : ''} |\n`;
  fs.writeFileSync(path.join(outDir, 'class-stats.md'), cmd);

  let und = '';
  for (const m of members.filter(dark).sort((a, b) => a.file.localeCompare(b.file)))
    und += `${m.file} | ${m.owner}.${m.name} | ${m.sig}\n`;
  fs.writeFileSync(path.join(outDir, 'dark-members.txt'), und);
  console.log(`tables written to ${outDir}`);
}

if (updateBaseline) {
  fs.writeFileSync(baselineFile, JSON.stringify(metrics, null, 2) + '\n');
  console.log(`baseline written to ${baselineFile}`);
}
else if (check) {
  const baseline = JSON.parse(fs.readFileSync(baselineFile, 'utf8'));
  const worse = Object.keys(baseline).filter((k) => metrics[k] > baseline[k]);
  const better = Object.keys(baseline).filter((k) => metrics[k] < baseline[k]);
  for (const k of better) console.log(`${k}: ${baseline[k]} -> ${metrics[k]} (better; run --update-baseline to keep it)`);
  for (const k of worse) console.error(`${k}: ${baseline[k]} -> ${metrics[k]} (worse)`);
  if (worse.length) {
    console.error('inventory: the JS API surface regressed; document/type/cover the new members, or update the baseline deliberately.');
    process.exit(1);
  }
  console.log('inventory: no regression against the baseline');
}
