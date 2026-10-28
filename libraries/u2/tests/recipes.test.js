import {test} from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, mkdirSync, readFileSync, readdirSync, writeFileSync} from 'node:fs';
import {createRequire, register} from 'node:module';
import {dirname, join, relative} from 'node:path';
import {fileURLToPath, pathToFileURL} from 'node:url';
import './dom-shim.js';

register('./dg-stub.mjs', import.meta.url);
const ts = createRequire(import.meta.url)('typescript');
const {Usage} = await import('datagrok-api/u2core');
const {expectOperable, mount, u2test} = await import('../src/testing/index.js');
const {ViewHost} = await import('../src/dg/index.js');

const root = fileURLToPath(new URL('..', import.meta.url));
const SAMPLES = join(root, 'docs', 'samples');
const OUT = join(root, '.test-build', 'docs');
const ENTRIES = {'@datagrok-libraries/u2': 'src/index', '@datagrok-libraries/u2/dg': 'src/dg/index',
  '@datagrok-libraries/u2/filter': 'src/core/filter/index', '@datagrok-libraries/u2/dg/filter': 'src/dg/filter/index',
  '@datagrok-libraries/u2/testing': 'src/testing/index'};
const NAMESPACES = {grok: 'datagrok-api/grok', ui: 'datagrok-api/ui', DG: 'datagrok-api/dg'};
const UNUSED = new Set([6133, 6138, 6192, 6196, 6198]);

const HEADLESS = {
  'docs/recipes/entity-browser.md › The shape': 'pages through grok.dapi.reports: a server',
  'docs/recipes/forms-and-details.md › The shape (metadata-driven, the preferred path)':
    'moleculeInput bridges the platform\'s molecule editor (ui.input.molecule)',
  'docs/recipes/shell-integration.md › An app that lives in a platform TableView':
    'needs a platform TableView: its grid, ribbon, dock manager and the shell\'s events',
  'docs/recipes/shell-integration.md › Anti-pattern': 'the anti-pattern: shown, never run',
  'docs/recipes/specs.md › Companion libraries register through a plugin':
    'ems runs against u2\'s built dist, not these sources',
  'docs/samples/table-view-chrome.ts': 'opens a platform TableView (grok.shell.addTableView)',
  'docs/samples/dataframe-binding.ts': 'docks into a platform TableView over a workspace table',
  'docs/samples/property-form.ts': 'builds DG.Property.fromOptions, which the node platform doubles do not have',
  '*.old.ts': 'the same job with the platform\'s ui.*, for the comparison: it needs the platform',
};

const VIEW = '(() => { const view = %s; __ctx.cleanup.push(() => ViewHost.of(view).dispose()); return view.root; })()';

const BLOCKS = {
  'README.md › Install': {run: true},
  'README.md › Testing': {kit: true,
    given: 'const test = (name: string, body: () => unknown) => void __ctx.tests.push(body);'},
  'docs/recipes/catalog-table.md › Composition': {run: 'host', given: `
    interface Campaign {
      code: string; friendlyName?: string; createDate: string; rowCount: number; status: string; modified: number;
      campaignFields?: Record<string, unknown>; template: {key: string, name: string};
    }
    const loadCampaigns = async (): Promise<Campaign[]> => [{code: 'KHL-12', friendlyName: 'PI3Kδ quick wins',
      createDate: '2026/10/19', rowCount: 120, status: 'In Progress', modified: 2, campaignFields: {Target: 'PI3K'},
      template: {key: 'khl', name: 'Kinase hit-to-lead'}}];
    const fieldsOf = (campaigns: Campaign[] | undefined) =>
      [...new Set((campaigns ?? []).flatMap((c) => Object.keys(c.campaignFields ?? {})))];
    const templatesWith = (field: string) => \`Kinase hit-to-lead (\${field})\`;
    const variantOf = (status: string): BadgeVariant => status === 'In Progress' ? 'accent' : 'default';
    const share = (_: Campaign): void => {};
    const remove = (_: Campaign): void => {};
    const editRights = (_: Campaign) => 'edit' as const;
    const open = (_: Campaign): void => {};
    const host = div();`},
  'docs/recipes/catalog-table.md › Keeping the choices where an app already keeps them': {run: true, given: `
    const readColumns = (text: string | null): string[] => text ? JSON.parse(text) : [];
    const fromLegacySort = (text: string | null) => text ?? '';`},
  'docs/recipes/entity-browser.md › The shape': {given: 'declare function reportRow(r: DG.UserReport): HTMLElement;'},
  'docs/recipes/forms-and-details.md › The shape (metadata-driven, the preferred path)': {given: `
    const props = [DG.Property.create('name', DG.TYPE.STRING, (x: any) => x.name, (x: any, v: unknown) => x.name = v),
      DG.Property.create('smiles', DG.TYPE.STRING, (x: any) => x.smiles, (x: any, v: unknown) => x.smiles = v)];
    const entity = {name: 'Chemists', smiles: 'CCO'} as unknown as DG.Group;
    const markDirty = (_: string): void => {};`},
  'docs/recipes/forms-and-details.md › Reading and writing a value': {run: 'personEditor()'},
  'docs/recipes/forms-and-details.md › Rows, picks and confirmations': {run: 'new Form().addAll([stages, order])',
    given: `
    const template = {stages: ['Design', 'Make', 'Test']};
    const functionNames = ['Descriptors', 'Solubility', 'Docking'];
    const saved = {order: ['Solubility', 'Descriptors'], hidden: ['Docking']};
    const deleteCampaign = async (_: string): Promise<void> => {};`},
  'docs/recipes/forms-and-details.md › Load, then build rows': {run: 'fieldsEditor()',
    given: 'const loadTemplate = async (name: string) => ({name, fields: [\'Compound\', \'IC50\']});'},
  'docs/recipes/forms-and-details.md › Dates and the text of a value': {run: 'dates'},
  'docs/recipes/forms-and-details.md › A docked editor panel': {run: 'panel', given: `
    const symbol = new TextInput({label: 'Symbol', value: 'Aib', required: true});
    const name = new TextInput({label: 'Name', value: 'α-Aminoisobutyric acid'});
    const polymer = new ChoiceInput({label: 'Polymer', items: ['PEPTIDE', 'RNA'], value: 'PEPTIDE'});
    const type = new ChoiceInput({label: 'Type', items: ['Backbone', 'Branch'], value: 'Backbone'});
    const id = new TextInput({label: 'ID', value: 'Aib-1'});
    const [line, background, text] = ['Line', 'Background', 'Text'].map((label) => new ColorInput({label, value: '#2083D5'}));
    const problem = signal<string | null>(null);
    const detailsProblem = signal<string | null>(null);
    const summary = 'Colors';
    const library = signal('Polytool');
    const onSave = (): void => {};
    const header = h3('Aib');
    const preview = div(['α-Aminoisobutyric acid']);
    const panel = div();`},
  'docs/recipes/list-item-rendering.md › Actions: one list, two surfaces': {run: 'list', given: `
    interface Report { id: string; description: string; jiraTicket?: string; createdOn: Date; }
    const copy = (_: string): void => {};
    const ticketUrl = (r: Report) => \`https://example.com/browse/\${r.jiraTicket}\`;
    const primaryLine = (r: Report) => span(r.description);`},
  'docs/recipes/list-item-rendering.md › Adaptive rendering: `u2-adaptive` + `u2-p*`': {run: 'list', given: `
    const list = new VirtualList<{assignee: string}>({itemHeight: 28, render: (x) => span(x.assignee)});
    const r = {assignee: 'ada'};`},
  'docs/recipes/shell-integration.md › The shape': {run: VIEW.replace('%s', 'reportsBrowserApp()'), given: `
    const buildReportsBrowser = () => ({content: div(['Reports']),
      ribbon: [[new TextInput({search: true, inline: true, placeholder: 'Search reports'})]], status: signal('3 reports')});`},
  'docs/recipes/shell-integration.md › An app that lives in a platform TableView': {given: `
    interface Session {
      recalc(cell: DG.GridCell): void; importFile(args: unknown): void; refresh(): void; editPermissions(): void;
      submit(): void; backToLanding(): void; apply(template: string): void;
      newChanges: ReadonlySignal<string | null>; canEdit: ReadonlySignal<boolean>; statusText: ReadonlySignal<string>;
    }
    declare const tv: DG.TableView;
    declare const landing: DG.ViewBase;
    declare const trackerDialog: Dialog;
    declare function campaignSession(tv: DG.TableView): Session;
    declare function detailsPanel(): Control;
    declare function loadTemplate(): Promise<string>;`},
  'docs/recipes/shell-integration.md › Routing': {run: VIEW.replace('%s', 'u2DemoApp(\'forms\')'), given: `
    const buildDemo = (o: {initial?: string}) =>
      ({content: div([o.initial ?? 'home']), shell: {current: signal('home')}, status: signal('Ready')});
    const leafForPath = (path?: string) => path === undefined ? undefined : {id: path};
    const bar = div();
    const demoRibbon = (_: {current: ReadonlySignal<string>}) => [div()];
    const pathOf = (leaf: string) => \`/\${leaf}\`;`},
  'docs/recipes/shell-integration.md › Anti-pattern': {given: `
    declare const search: TextInput;
    declare const split: HTMLElement;
    declare const statusStrip: HTMLElement;
    declare function refresh(): void;`},
  'docs/recipes/specs.md › Companion libraries register through a plugin': {},
};

const SAMPLE_RUNS = {
  'app-view': (m, ctx) => viewRoot(m.projectsApp(), ctx),
  'app-home': (m, ctx) => viewRoot(m.myApp(), ctx),
};

function viewRoot(view, ctx) {
  ctx.cleanup.push(() => ViewHost.of(view).dispose());
  return view.root;
}

const rel = (file) => relative(root, file).split('\\').join('/');
const read = (file) => readFileSync(file, 'utf8').replace(/\r\n/g, '\n');

function markdownBlocks() {
  const docs = [join(root, 'README.md'), ...readdirSync(join(root, 'docs', 'recipes')).filter((f) => f.endsWith('.md'))
    .map((f) => join(root, 'docs', 'recipes', f))];
  const blocks = [];
  for (const file of docs) {
    const lines = read(file).split('\n');
    const counts = new Map();
    let heading = '';
    for (let i = 0; i < lines.length; i++) {
      const h = /^#{1,6} (.*)$/.exec(lines[i]);
      if (h !== null)
        heading = h[1];
      const fence = /^```(\w*)\s*$/.exec(lines[i]);
      if (fence === null)
        continue;
      const end = lines.indexOf('```', i + 1);
      if (fence[1] === 'ts' || fence[1] === 'typescript') {
        const n = (counts.get(heading) ?? 0) + 1;
        counts.set(heading, n);
        blocks.push({key: `${rel(file)} › ${heading}${n > 1 ? ` #${n}` : ''}`, file: rel(file), line: i + 2,
          code: lines.slice(i + 1, end).join('\n')});
      }
      i = end;
    }
  }
  return blocks;
}

function exportsOf(entry, seen = new Set()) {
  const file = join(root, `${entry}.ts`);
  if (seen.has(file))
    return [];
  seen.add(file);
  const sf = ts.createSourceFile(file, read(file), ts.ScriptTarget.Latest, true);
  const names = [];
  for (const st of sf.statements) {
    if (ts.isExportDeclaration(st)) {
      if (st.exportClause !== undefined && ts.isNamedExports(st.exportClause))
        names.push(...st.exportClause.elements.map((e) => e.name.text));
      else if (st.exportClause === undefined && st.moduleSpecifier !== undefined)
        names.push(...exportsOf(join(dirname(entry), st.moduleSpecifier.text.replace(/\.js$/, '')), seen));
    }
    else if (st.modifiers?.some((m) => m.kind === ts.SyntaxKind.ExportKeyword)) {
      if (ts.isVariableStatement(st))
        names.push(...st.declarationList.declarations.flatMap((d) => bound(d.name)));
      else if (st.name !== undefined)
        names.push(st.name.text);
    }
  }
  return names;
}

function bound(name) {
  return ts.isIdentifier(name) ? [name.text] : name.elements.flatMap((e) => ts.isOmittedExpression(e) ? [] : bound(e.name));
}

const PROVIDER = new Map();
for (const [spec, entry] of Object.entries(ENTRIES).slice(0, 2)) {
  for (const name of exportsOf(entry))
    if (!PROVIDER.has(name))
      PROVIDER.set(name, spec);
}

function declaredAt(sf) {
  const out = new Set();
  for (const st of sf.statements) {
    if (ts.isImportDeclaration(st)) {
      const clause = st.importClause;
      if (clause?.name)
        out.add(clause.name.text);
      if (clause?.namedBindings && ts.isNamespaceImport(clause.namedBindings))
        out.add(clause.namedBindings.name.text);
      else if (clause?.namedBindings)
        clause.namedBindings.elements.forEach((e) => out.add(e.name.text));
    }
    else if (ts.isVariableStatement(st))
      st.declarationList.declarations.flatMap((d) => bound(d.name)).forEach((n) => out.add(n));
    else if (st.name !== undefined && ts.isIdentifier(st.name))
      out.add(st.name.text);
  }
  return out;
}

function referenced(sf) {
  const out = new Set();
  const skip = (node) => (ts.isPropertyAccessExpression(node.parent) && node.parent.name === node) ||
    (ts.isQualifiedName(node.parent) && node.parent.right === node) ||
    ((ts.isPropertyAssignment(node.parent) || ts.isPropertySignature(node.parent) ||
      ts.isMethodDeclaration(node.parent) || ts.isPropertyDeclaration(node.parent)) && node.parent.name === node);
  const visit = (node) => {
    if (ts.isIdentifier(node) && !skip(node))
      out.add(node.text);
    ts.forEachChild(node, visit);
  };
  visit(sf);
  return out;
}

function blockSources(block, spec) {
  const given = (spec.given ?? '').replace(/^\n/, '').split('\n').map((l) => l.replace(/^ {4}/, '')).join('\n');
  const run = typeof spec.run === 'string' ? spec.run : 'undefined';
  const sf = ts.createSourceFile('block.ts', `${given}\n${block.code}\n${run};`, ts.ScriptTarget.Latest, true);
  const own = declaredAt(sf);
  const auto = new Map();
  for (const name of referenced(sf)) {
    if (own.has(name))
      continue;
    const from = NAMESPACES[name] ? null : PROVIDER.get(name);
    if (from)
      auto.set(from, [...auto.get(from) ?? [], name]);
  }
  const head = [...auto].map(([from, names]) => `import {${names.join(', ')}} from '${from}';`);
  for (const [name, from] of Object.entries(NAMESPACES))
    if (!own.has(name) && referenced(sf).has(name))
      head.push(`import * as ${name} from '${from}';`);
  head.push('declare const __ctx: {tests: (() => unknown)[], cleanup: (() => void)[]};');
  const code = ts.createSourceFile('code.ts', block.code, ts.ScriptTarget.Latest, true);
  const imports = code.statements.filter((s) => ts.isImportDeclaration(s)).map((s) => s.getText(code));
  const body = code.statements.filter((s) => !ts.isImportDeclaration(s)).map((s) =>
    s.getFullText(code).slice(0, s.getStart(code) - s.getFullStart()) + s.getText(code).replace(/^export\s+/, ''))
    .join('');
  const compiled = [...head, given].join('\n');
  return {
    compile: `${compiled}\n${block.code}\n`,
    offset: compiled.split('\n').length,
    run: [...imports, ...head.slice(0, -1), 'export default function __run(__ctx) {', given, body,
      `return (${run});`, '}'].join('\n'),
  };
}

function programOver(virtual) {
  const config = join(SAMPLES, 'tsconfig.json');
  const parsed = ts.getParsedCommandLineOfConfigFile(config, {}, {...ts.sys, onUnRecoverableConfigFileDiagnostic: (d) => {
    throw new Error(ts.flattenDiagnosticMessageText(d.messageText, '\n'));
  }});
  const options = {...parsed.options, paths: {...parsed.options.paths,
    '@datagrok-libraries/ems/src/dg': ['../../../ems/dist/src/dg/index.d.ts']}};
  const host = ts.createCompilerHost(options, true);
  const getSourceFile = host.getSourceFile.bind(host);
  host.getSourceFile = (name, version, onError, create) => virtual.has(name) ?
    ts.createSourceFile(name, virtual.get(name), version, true) : getSourceFile(name, version, onError, create);
  const fileExists = host.fileExists.bind(host);
  host.fileExists = (name) => virtual.has(name) || fileExists(name);
  const readFile = host.readFile.bind(host);
  host.readFile = (name) => virtual.get(name) ?? readFile(name);
  return {program: ts.createProgram([...parsed.fileNames, ...virtual.keys()], options, host), samples: parsed.fileNames};
}

function messages(program, sf, keep) {
  return [...program.getSyntacticDiagnostics(sf), ...program.getSemanticDiagnostics(sf)].filter(keep).map((d) => {
    const at = d.start === undefined ? 0 : sf.getLineAndCharacterOfPosition(d.start).line;
    return {line: at + 1, text: ts.flattenDiagnosticMessageText(d.messageText, ' ')};
  });
}

const slug = (key) => key.replace(/[^a-zA-Z0-9]+/g, '-').replace(/^-|-$/g, '');

function emit(name, source) {
  const js = ts.transpileModule(source, {compilerOptions: {module: ts.ModuleKind.ESNext,
    target: ts.ScriptTarget.ES2022, useDefineForClassFields: false}}).outputText
    .replace(/from (['"])(@datagrok-libraries\/u2[\w/]*)\1/g, (m, q, spec) =>
      ENTRIES[spec] ? `from '${pathToFileURL(join(root, '.test-build', `${ENTRIES[spec]}.js`)).href}'` : m);
  mkdirSync(OUT, {recursive: true});
  const file = join(OUT, `${name}.mjs`);
  writeFileSync(file, js);
  return pathToFileURL(file).href;
}

async function strict(body) {
  const mode = Usage.mode;
  const said = [];
  const {warn, error} = console;
  const hear = (log) => (...args) => {
    if (typeof args[0] === 'string' && args[0].startsWith('u2'))
      said.push(args.map(String).join(' '));
    log(...args);
  };
  Usage.mode = 'throw';
  console.warn = hear(warn);
  console.error = hear(error);
  try {
    await body();
  }
  finally {
    Usage.mode = mode;
    Object.assign(console, {warn, error});
  }
  assert.deepEqual(said, [], 'u2 said nothing to the console');
}

async function live(build) {
  const ctx = {tests: [], cleanup: []};
  await strict(() => u2test(() => {
    try {
      const built = build(ctx);
      if (built !== undefined)
        expectOperable(mount(() => built));
    }
    finally {
      for (const done of ctx.cleanup.splice(0))
        done();
    }
  })());
}

const blocks = markdownBlocks();
const virtual = new Map();
const sources = new Map();
for (const block of blocks) {
  const spec = BLOCKS[block.key];
  if (spec === undefined)
    continue;
  const s = blockSources(block, spec);
  const file = join(root, 'docs', '.blocks', `${slug(block.key)}.ts`);
  sources.set(block.key, {...s, file});
  virtual.set(file, s.compile);
}
const {program, samples} = programOver(virtual);
const sampleIds = samples.map((f) => rel(f));
const headless = (key) => HEADLESS[key] ?? (key.endsWith('.old.ts') ? HEADLESS['*.old.ts'] : undefined);

test('every TypeScript block of the README and the recipes is listed, with what it runs or why it cannot', () => {
  assert.deepEqual(blocks.map((b) => b.key).filter((k) => BLOCKS[k] === undefined), []);
  assert.deepEqual(Object.keys(BLOCKS).filter((k) => !blocks.some((b) => b.key === k)), [], 'no entry for a gone block');
  assert.deepEqual(Object.keys(HEADLESS).filter((k) => k !== '*.old.ts' && BLOCKS[k] === undefined &&
    !sampleIds.includes(k)), [], 'no headless mark for a gone block or sample');
  for (const [key, spec] of Object.entries(BLOCKS))
    assert.ok((spec.run !== undefined || spec.kit) !== (HEADLESS[key] !== undefined), `${key}: runs, or says why not`);
});

test('every TypeScript block of the README and the recipes compiles against u2\'s sources', () => {
  const broken = [];
  for (const block of blocks.filter((b) => sources.has(b.key))) {
    const {file, offset} = sources.get(block.key);
    for (const m of messages(program, program.getSourceFile(file), (d) => !UNUSED.has(d.code))) {
      const where = m.line > offset ? `${block.file}:${block.line + m.line - offset - 1}` : `${block.key} (given)`;
      broken.push(`${where}: ${m.text}`);
    }
  }
  assert.deepEqual(broken, []);
});

test('every sample compiles against u2\'s sources, unused names included', () => {
  assert.ok(samples.length >= 38, `${samples.length} samples`);
  const broken = samples.flatMap((f) => messages(program, program.getSourceFile(f), () => true)
    .map((m) => `${rel(f)}:${m.line}: ${m.text}`));
  assert.deepEqual(broken, []);
});

test('every sample has its board text: front matter, the live demo it says, its CSS and its comparison', () => {
  const wrong = [];
  for (const file of samples.filter((f) => !f.endsWith('.old.ts'))) {
    const id = rel(file).slice('docs/samples/'.length, -3);
    const md = join(SAMPLES, `${id}.md`);
    if (!existsSync(md)) {
      wrong.push(`${id}: no ${id}.md`);
      continue;
    }
    const front = /^---\n([\s\S]*?)\n---\n/.exec(read(md));
    const meta = Object.fromEntries((front?.[1] ?? '').split('\n').map((l) => /^(\w+):\s*(.*)$/.exec(l))
      .filter((m) => m !== null).map((m) => [m[1], m[2]]));
    const code = read(file);
    for (const key of ['title', 'group', 'summary', 'live'])
      if (!meta[key])
        wrong.push(`${id}.md: no ${key}`);
    if (!read(md).includes('\n## What it takes\n'))
      wrong.push(`${id}.md: no "What it takes"`);
    if ((meta.live === 'yes') !== /^export function render\(host: HTMLElement\)/m.test(code))
      wrong.push(`${id}: live is ${meta.live}, render(host) ${code.includes('function render(') ? 'is' : 'is not'} there`);
    if (meta.live === 'yes' && !(code.indexOf('\n// Live demo\n') >= 0 && code.indexOf('// Live demo') < code.indexOf('function render(')))
      wrong.push(`${id}: the demo goes under a "// Live demo" line`);
    if (existsSync(join(SAMPLES, `${id}.css`)) !== code.includes(`import './${id}.css';`))
      wrong.push(`${id}: ${id}.css and its import go together`);
    if (existsSync(join(SAMPLES, `${id}.old.ts`)) !== (meta.old !== undefined))
      wrong.push(`${id}: ${id}.old.ts and "old:" go together`);
  }
  assert.deepEqual(wrong, []);
});

for (const block of blocks) {
  const spec = BLOCKS[block.key];
  if (spec === undefined || headless(block.key) !== undefined)
    continue;
  test(`recipe runs headless, strict, leaving nothing live: ${block.key}`, async () => {
    const mod = await import(emit(slug(block.key), sources.get(block.key).run));
    if (spec.kit) {
      const ctx = {tests: [], cleanup: []};
      await strict(async () => {
        mod.default(ctx);
        assert.ok(ctx.tests.length > 0, 'the block declares its tests');
        for (const body of ctx.tests)
          await body();
      });
      return;
    }
    await live((ctx) => mod.default(ctx));
  });
}

for (const file of samples) {
  const key = rel(file);
  const id = key.slice('docs/samples/'.length, -3);
  if (headless(key) !== undefined)
    continue;
  test(`sample runs headless, strict, leaving nothing live, the probe clean: ${id}`, async () => {
    const mod = await import(emit(`sample-${id}`, read(file).replace(/^import '\.\/[\w-]+\.css';\n/m, '')));
    await live((ctx) => {
      if (SAMPLE_RUNS[id] !== undefined)
        return SAMPLE_RUNS[id](mod, ctx);
      assert.equal(typeof mod.render, 'function', `${id} is live: it exports render(host)`);
      const host = document.createElement('div');
      mod.render(host);
      return host;
    });
  });
}
