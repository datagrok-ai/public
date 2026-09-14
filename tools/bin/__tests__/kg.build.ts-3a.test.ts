/// `grok kg build` WO-3a (build-plan.md): the annotation parser, and the ts-packages / ts-functions
/// extractors against the mini monorepo under fixtures/kg/build, whose type files copy the real ones.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {parseHeaderLines, parseParam, parseFunctionHeaders, parseScriptHeader, parseQueryHeaders} from '../utils/kg/build/annotations';
import {Emitter} from '../utils/kg/build/emitter';
import {loadTypeSystem} from '../utils/kg/types';
import {currentDir} from '../utils/kg/build/write';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const realKg = path.resolve(fixture, '..', '..', '..', '..', '..', '..', '..', 'core', 'docs', 'knowledge-graph');

async function build(): Promise<{manifest: any, rows: (file: string) => any[], problems: Record<string, string[]>}> {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-3a-'));
  fs.cpSync(fixture, repo, {recursive: true});
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'ts-packages,ts-functions', db: false, output: 'json'});
    expect(error.mock.calls).toEqual([]);
    expect(process.exitCode).toBeUndefined();
    const out = currentDir(path.join(repo, '.kg'))!;
    const rows = (file: string) => {
      const p = path.join(out, file.startsWith('reports/') ? file : `data/${file}.jsonl`);
      return fs.existsSync(p) ? fs.readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l)) : [];
    };
    return {manifest: JSON.parse(String(log.mock.calls[0][0])), rows, problems: JSON.parse(fs.readFileSync(path.join(out, 'reports', 'problems.json'), 'utf8'))};
  } finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

const graph = build();
const byId = (rows: any[], id: string) => rows.find((r) => r.id === id);

describe('annotation parser (build-plan.md WO-3a)', () => {
  it('keeps every key verbatim: meta.* into meta, the rest under keys, repeated keys as lists', () => {
    const h = parseHeaderLines([
      '//name: Recalculate Coordinates [x]', '//description: Recalculates 2D coordinates', '//help-url: /help/chem', '//top-menu: Chem | Transform',
      '//condition: true', '//editor-for: AddNewColumn', '//connection: Chembl', '//environment: chem', '//test: expect(1, 1)', '//test: expect(2, 2)',
      '//sample: a.py', '//reference: https://x', '//friendlyName: Recalc', '//feature: ~visualize/x', '//tags: app, panel', '//meta.role: transform',
      '//meta.cache: all', '//meta.cache.invalidateOn: 0 0 * * *', '//meta.method_info.year: 2024', '// TODO: unknown key', '//// not a header',
    ], 7)!;
    expect(h.line).toBe(7);
    expect(h.name).toBe('Recalculate Coordinates');
    expect(h.description).toBe('Recalculates 2D coordinates');
    expect(h.tags).toEqual(['app', 'panel']);
    expect(h.meta).toEqual({role: 'transform', cache: 'all', 'cache.invalidateOn': '0 0 * * *', 'method_info.year': '2024'});
    expect(h.keys).toMatchObject({'help-url': ['/help/chem'], 'top-menu': ['Chem | Transform'], condition: ['true'], 'editor-for': ['AddNewColumn'],
      connection: ['Chembl'], environment: ['chem'], test: ['expect(1, 1)', 'expect(2, 2)'], sample: ['a.py'], reference: ['https://x'],
      friendlyName: ['Recalc'], feature: ['~visualize/x'], 'meta.role': ['transform'], TODO: ['unknown key']});
    expect(parseHeaderLines(['// TODO: nothing known', '// https://example.com'], 1)).toBeNull();
    expect(parseHeaderLines(['//meta.only: x'], 1)).not.toBeNull();
  });

  it('parses input and output options in braces, defaults and bracketed descriptions', () => {
    expect(parseParam('column molecules { semType: Molecule }')).toEqual({type: 'column', name: 'molecules', options: {semType: 'Molecule'}});
    expect(parseParam('string method { choices: ["OCL","CoordGen"] }')).toEqual({type: 'string', name: 'method', options: {choices: '["OCL","CoordGen"]'}});
    expect(parseParam('string result { semType: Macromolecule; units: helm }')).toMatchObject({options: {semType: 'Macromolecule', units: 'helm'}});
    expect(parseParam('string path {meta.url: true; optional: true}')).toMatchObject({options: {'meta.url': 'true', optional: 'true'}});
    expect(parseParam('double pH = 7.4 {caption: pH}')).toEqual({type: 'double', name: 'pH', default: '7.4', options: {caption: 'pH'}});
    expect(parseParam('string target = "CHEMBL1827" [ChEMBL target identifier, e.g. CHEMBL1827]')).toEqual({type: 'string', name: 'target', default: '"CHEMBL1827"', options: {}, description: 'ChEMBL target identifier, e.g. CHEMBL1827'});
    expect(parseParam('dataframe result {action:join(table)}')).toMatchObject({options: {action: 'join(table)'}});
    expect(parseParam('list x = [1, 2]')).toMatchObject({default: '[1, 2]', options: {}});
    expect(parseParam('object module')).toEqual({type: 'object', name: 'module', options: {}});
  });

  it('accepts function, arrow and class-method starts, so detectors.js parses, and hands each block its body', () => {
    const blocks = parseFunctionHeaders([
      'class ChemPackageDetectors extends DG.Package {', '  static likelyNames = [];', '',
      '  //name: detectMolecules', '  //meta.role: semTypeDetector', '  //input: column col', '  detectMolecules(col) {', '    return DG.SEMTYPE.MOLECULE;', '  }', '',
      '  //meta.role: semTypeDetector', '  static async detectReactions(col: DG.Column): Promise<string | null> {', '    return null;', '  }', '}',
      '//name: Scaffold Tree', '//meta.role: viewer', 'export function scaffoldTreeViewer() : any {', '}',
      '//name: Arrow', 'export const arrow = async (x: number): Promise<void> => {', '};',
      '// not a header', 'if (x) {', '}',
      '//name: gap', '', 'export function gap() {}',
    ].join('\n'));
    expect(blocks.map((b) => [b.header.name ?? b.declaration.name, b.declaration.kind, b.declaration.owner, b.declaration.line])).toEqual([
      ['detectMolecules', 'method', 'ChemPackageDetectors', 7], ['detectReactions', 'method', 'ChemPackageDetectors', 12],
      ['Scaffold Tree', 'function', undefined, 18], ['Arrow', 'arrow', undefined, 21],
    ]);
    expect(blocks[0].body).toContain('return DG.SEMTYPE.MOLECULE;');
    expect(blocks[0].body).not.toContain('detectReactions');
    expect(blocks[0].header.line).toBe(4);
  });

  it('reads the # prefix of scripts and the -- prefix of queries, several queries per file', () => {
    const script = parseScriptHeader('\n#name: Calculate logD\n#language: python\n#input: column molecules {semType: Molecule}\n## a comment\n#output: dataframe result\n\nimport pandas\n#name: not a header\n', 'python')!;
    expect(script).toMatchObject({name: 'Calculate logD', line: 2, keys: {language: ['python']}});
    expect(script.inputs).toEqual([{type: 'column', name: 'molecules', options: {semType: 'Molecule'}}]);
    expect(script.outputs).toEqual([{type: 'dataframe', name: 'result', options: {}}]);
    expect(parseScriptHeader('import os\n#name: late\n', 'python')).toBeNull();
    const queries = parseQueryHeaders('--name: a\n--connection: C\nselect 1\n--end\n\n-- name: b\n--input: int n = 1\nselect 2\n--comment: only\nselect 3\n');
    expect(queries.map((q) => [q.name, q.line, q.keys.connection])).toEqual([['a', 1, ['C']], ['b', 6, undefined]]);
  });

  it('parses a decorator-lowered package.g.ts sample the way check.ts does, name verbatim', () => {
    const blocks = parseFunctionHeaders([
      '//name: Recalculate Coordinates', '//description: Recalculates 2D coordinates for molecules', '//input: dataframe table ', '//input: column molecules { semType: Molecule }',
      '//input: string method { choices: ["OCL","CoordGen"] }', '//input: bool join = true ', '//output: column result', '//meta.role: transform',
      '//top-menu: Chem | Transform | Recalculate Coordinates...',
      'export async function recalculateCoords(table: DG.DataFrame, molecules: DG.Column, method: string, join: boolean) : Promise<any> {',
      '  return await PackageFunctions.recalculateCoords(table, molecules, method, join);', '}',
    ].join('\n'));
    expect(blocks).toHaveLength(1);
    const {header, declaration} = blocks[0];
    expect(header.name).toBe('Recalculate Coordinates');
    expect(declaration).toEqual({line: 10, name: 'recalculateCoords', kind: 'function'});
    expect(header.inputs.map((p) => [p.type, p.name, p.default, p.options])).toEqual([
      ['dataframe', 'table', undefined, {}], ['column', 'molecules', undefined, {semType: 'Molecule'}],
      ['string', 'method', undefined, {choices: '["OCL","CoordGen"]'}], ['bool', 'join', 'true', {}],
    ]);
    expect(header.meta).toEqual({role: 'transform'});
    expect(header.keys['top-menu']).toEqual(['Chem | Transform | Recalculate Coordinates...']);
  });
});

describe('ts-packages extractor (build-plan.md WO-3a)', () => {
  it('copies the fixture type files from the real ones, or skips when the monorepo is not around', () => {
    if (!fs.existsSync(path.join(realKg, 'schema.yaml'))) return;
    const files = (root: string) => ['schema.yaml', ...['nodes', 'edges'].flatMap((d) => fs.readdirSync(path.join(root, d)).map((f) => `${d}/${f}`))].sort();
    const fixtureKg = path.join(fixture, KG_DIR);
    expect(files(fixtureKg)).toEqual(files(realKg));
    for (const f of files(fixtureKg)) expect(fs.readFileSync(path.join(fixtureKg, f), 'utf8'), f).toBe(fs.readFileSync(path.join(realKg, f), 'utf8'));
  });

  it('reads package.json into package and library nodes with the type-file members', async () => {
    const {rows, manifest} = await graph;
    expect(manifest.sources).toEqual({'ts-functions': 'partial(2 rejected)', 'ts-packages': 'ok'});
    expect(byId(rows('nodes/package'), 'pkg:Demo')).toEqual({
      id: 'pkg:Demo', type: 'package', name: 'Demo', author: 'Jane Dev', batch: expect.any(String), category: 'Cheminformatics', description: 'The fixture package with every folder.',
      friendly_name: 'Demo', language: 'ts', npm: '@datagrok/demo', path: 'public/packages/Demo', provenance: 'registry', service: true, settings: ['Sketcher', 'TemplatesPath'],
      source_layer: 'public', sources: ['common/openchemlib-full.js'], status: 'active', version: '1.2.3', visibility: 'public',
    });
    expect(byId(rows('nodes/package'), 'pkg:Plain')).toMatchObject({language: 'ts', service: false, npm: '@datagrok/plain'});
    expect(rows('nodes/library')).toEqual([
      expect.objectContaining({id: 'lib:js-api', type: 'library', name: 'js-api', npm: 'datagrok-api', version: '1.27.11', language: 'ts', path: 'public/js-api', provenance: 'registry'}),
      expect.objectContaining({id: 'lib:utils', name: 'utils', npm: '@datagrok-libraries/utils', version: '4.7.9', path: 'public/libraries/utils'}),
    ]);
  });

  it('classifies dependencies into depends-on with kind and range; other npm names and unknown targets are left out', async () => {
    const {rows, problems} = await graph;
    expect(rows('edges/depends-on').map((e) => [e.from, e.to, e.kind, e.range])).toEqual([
      ['lib:utils', 'lib:js-api', 'runtime', '^1.27.7'],
      ['pkg:Demo', 'lib:js-api', 'runtime', '^1.27.7'],
      ['pkg:Demo', 'lib:utils', 'optional', '^4.7.9'],
      ['pkg:Demo', 'lib:utils', 'runtime', '^4.7.9'],
      ['pkg:Demo', 'pkg:Plain', 'dev', '^0.1.0'],
      ['pkg:Plain', 'lib:js-api', 'runtime', '^1.27.7'],
    ]);
    expect(rows('edges/depends-on')[1]).toMatchObject({derived_by: 'registry', confidence: 1, evidence: ['public/packages/Demo/package.json']});
    expect(problems.unresolved_ids).toContain('public/packages/Demo/package.json: peerDependencies @datagrok-libraries/nowhere is not under public/packages or public/libraries');
  });

  it('declares the semantic types of package.json with declared_in, and the package declares everything it owns', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/semantic-type'), 'semtype:DemoId')).toMatchObject({name: 'DemoId', description: 'A demo identifier', language: 'other', declared_in: 'pkg:Demo', provenance: 'registry'});
    expect(rows('edges/declared_in')).toEqual([expect.objectContaining({type: 'ref', name: 'declared_in', from: 'semtype:DemoId', to: 'pkg:Demo'})]);
    const owned = rows('edges/declares').filter((e) => e.from === 'pkg:Demo').map((e) => e.to);
    for (const id of ['semtype:DemoId', 'func:Demo:To HELM', 'func:Demo:detectMolecules', 'func:Demo:Calculate logD', 'func:Demo:users', 'conn:Demo:Demo',
      'env:Demo:demo-env', 'container:Demo:demo', 'file:public/packages/Demo/src/package.g.ts', 'file:public/packages/Demo/detectors.js'])
      expect(owned).toContain(id);
    expect(rows('edges/declares').find((e) => e.from === 'pkg:Demo' && e.to === 'func:Demo:To HELM')).toMatchObject({derived_by: 'registry', evidence: ['public/packages/Demo/src/package.g.ts']});
  });
});

describe('ts-functions extractor (build-plan.md WO-3a)', () => {
  it('picks the subtype by role precedence: widgets,panel is a panel; viewer,panel a viewer; adminApp,app an app; every role stays in roles', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/panel'), 'func:Demo:Molecule Panel')).toMatchObject({type: 'panel', roles: ['panel', 'widgets']});
    expect(byId(rows('nodes/viewer'), 'func:Demo:Demo Viewer')).toMatchObject({type: 'viewer', roles: ['panel', 'viewer'], trellisable: true, icon: 'files/icons/viewer.svg', builtin: false});
    expect(byId(rows('nodes/app'), 'func:Demo:Demo App')).toMatchObject({type: 'app', roles: ['adminApp', 'app'], admin: true, browse_path: 'Chem | Demo', icon: 'files/icons/demo.svg'});
    expect(byId(rows('nodes/app'), 'func:Plain:Plain App')).toMatchObject({admin: false, browse_path: 'Misc', tags: ['app'], roles: ['app']});
    expect(rows('nodes/lifecycle-hook').map((f) => [f.id, f.phase, f.immediate])).toEqual([['func:Demo:autostart', 'autostart', true], ['func:Demo:init', 'init', undefined]]);
  });

  it('fills the subtype properties from meta and drops them from the bag; a role missing its required member stays a function', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/cell-renderer'), 'func:Demo:Molecule Renderer')).toMatchObject({cell_type: 'Molecule', column_tags: ['quality=Molecule', 'foo=bar'], meta: {}});
    expect(byId(rows('nodes/file-handler'), 'func:Demo:Import SDF')).toMatchObject({extensions: ['sdf', 'mol'], direction: 'import', meta: {}});
    expect(byId(rows('nodes/file-viewer'), 'func:Demo:Preview MOL')).toMatchObject({extensions: ['mol', 'mol2'], check: 'Demo:checkMol', meta: {}});
    expect(byId(rows('nodes/panel'), 'func:Demo:Molecule Panel')).toMatchObject({condition: 'true', target_type: 'string', target_semtype: 'semtype:Molecule', meta: {}});
    expect(byId(rows('nodes/filter'), 'func:Demo:Substructure Filter')).toMatchObject({semtype: 'semtype:Molecule', primary: true, columnless: true, meta: {}});
    expect(byId(rows('nodes/editor'), 'func:Demo:Column Editor')).toMatchObject({meta: {'editor-for': 'AddNewColumn'}});
    expect(byId(rows('nodes/script-handler'), 'func:Demo:Demo Handler')).toMatchObject({script_language: 'demo', extensions: ['dm'], comment_start: '#', meta: {}});
    expect(byId(rows('nodes/sem-type-detector'), 'func:Demo:detectSequences')).toMatchObject({skip_test: true, meta: {}});
    expect(byId(rows('nodes/function'), 'func:Demo:Broken Renderer')).toMatchObject({type: 'function', roles: ['cellRenderer']});
    expect(rows('nodes/cell-renderer').map((f) => f.id)).toEqual(['func:Demo:Molecule Renderer']);
  });

  it('lifts the common members and keeps every other key in meta verbatim', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/function'), 'func:Demo:To HELM')).toEqual({
      id: 'func:Demo:To HELM', type: 'function', name: 'To HELM', batch: expect.any(String), cache: 'all', demo_path: 'Bioinformatics | To HELM',
      description: 'Converts a sequence to HELM', help_url: 'https://datagrok.ai/help/domains/bio', input_types: ['dataframe', 'column', 'string', 'int'], language: 'ts', line: 119,
      meta: {'cache.invalidateOn': '0 0 * * *', vectorFunc: 'true'}, output_types: ['column'], package: 'pkg:Demo', path: 'public/packages/Demo/src/package.g.ts',
      provenance: 'annotation', signature: '(dataframe table, column sequence: Macromolecule, string method, int n) -> column', source_layer: 'public', status: 'active',
      top_menu: 'Bio | Convert | To HELM...', visibility: 'public',
    });
    expect(byId(rows('nodes/function'), 'func:Demo:toMolfile')).toMatchObject({name: 'toMolfile', signature: '(string mol) -> string'});
    expect(rows('nodes/source-file').map((f) => [f.id, f.generated])).toEqual([
      ['file:public/packages/Demo/detectors.js', false], ['file:public/packages/Demo/package.js', false], ['file:public/packages/Demo/src/package.g.ts', true],
      ['file:public/packages/Demo/src/package.ts', false], ['file:public/packages/Demo/src/utils.ts', false],
      ['file:public/packages/Plain/src/package-test.ts', false], ['file:public/packages/Plain/src/package.js', false], ['file:public/packages/Plain/src/package.ts', false],
    ]);
  });

  it('falls back to package.ts and package-test.ts when there is no package.g.ts', async () => {
    const {rows} = await graph;
    expect(rows('edges/declares').filter((e) => e.from.startsWith('file:public/packages/Plain/')).map((e) => [e.from, e.to])).toEqual([
      ['file:public/packages/Plain/src/package-test.ts', 'func:Plain:test'],
      ['file:public/packages/Plain/src/package.js', 'func:Plain:srcJs'],
      ['file:public/packages/Plain/src/package.ts', 'func:Plain:Plain App'],
      ['file:public/packages/Plain/src/package.ts', 'func:Plain:helper'],
    ]);
  });

  it('turns //feature: into a rung-1 marker claim for the file and a direct is-implemented-in edge, stubbing the feature', async () => {
    const {rows} = await graph;
    expect(rows('reports/claims.jsonl')).toEqual([
      {feature: 'domains/bio', file: 'public/packages/Demo/scripts/run.js', line: 1, props: {}, rung: 1, source: 'marker'},
      {feature: 'domains/bio', file: 'public/packages/Demo/src/package.g.ts', line: 105, props: {}, rung: 1, source: 'marker'},
    ]);
    expect(rows('edges/is-implemented-in')).toEqual([
      expect.objectContaining({from: 'domains/bio', to: 'func:Demo:RunJs', derived_by: 'annotation', confidence: 1, evidence: ['public/packages/Demo/scripts/run.js']}),
      expect.objectContaining({from: 'domains/bio', to: 'func:Demo:To HELM', evidence: ['public/packages/Demo/src/package.g.ts']}),
    ]);
    expect(byId(rows('nodes/feature'), 'domains/bio')).toMatchObject({status: 'proposed', provenance: 'annotation'});
  });

  it('declares a function from its file and, for a method, from the class declaration stub', async () => {
    const {rows} = await graph;
    const decl = 'decl:public/packages/Demo/detectors.js#DemoPackageDetectors';
    expect(byId(rows('nodes/declaration'), decl)).toMatchObject({name: 'DemoPackageDetectors', kind: 'class', language: 'js', path: 'public/packages/Demo/detectors.js', status: 'proposed', provenance: 'ast'});
    expect(rows('edges/declares').filter((e) => e.to === 'func:Demo:detectMolecules').map((e) => [e.from, e.derived_by])).toEqual([
      [decl, 'ast'], ['file:public/packages/Demo/detectors.js', 'ast'], ['pkg:Demo', 'registry'],
    ]);
  });

  it('emits targets-semtype with the four roles, detection from meta.semType or from what the body returns or assigns', async () => {
    const {rows} = await graph;
    expect(rows('edges/targets-semtype').map((e) => [e.from, e.to, e.role, e.derived_by, e.confidence])).toEqual([
      ['func:Demo:Calculate logD', 'semtype:Molecule', 'consumes', 'annotation', 1],
      ['func:Demo:Molecule Panel', 'semtype:Molecule', 'consumes', 'annotation', 1],
      ['func:Demo:Molecule Renderer', 'semtype:Molecule', 'renders', 'annotation', 1],
      ['func:Demo:Substructure Filter', 'semtype:Molecule', 'filters', 'annotation', 1],
      ['func:Demo:To HELM', 'semtype:Macromolecule', 'consumes', 'annotation', 1],
      ['func:Demo:To HELM', 'semtype:Macromolecule', 'produces', 'annotation', 1],
      ['func:Demo:detectCountries', 'semtype:demo-country', 'detects', 'ast', 0.9],
      ['func:Demo:detectFlags', 'semtype:flag', 'detects', 'ast', 0.9],
      ['func:Demo:detectImages', 'semtype:BinaryImage', 'detects', 'ast', 0.9],
      ['func:Demo:detectImages', 'semtype:Text', 'detects', 'ast', 0.9],
      ['func:Demo:detectMolecules', 'semtype:Molecule', 'detects', 'ast', 0.9],
      ['func:Demo:detectSequences', 'semtype:Macromolecule', 'detects', 'annotation', 1],
      ['func:Demo:toMolfile', 'semtype:Molecule', 'produces', 'annotation', 1],
    ]);
    expect(rows('nodes/semantic-type').map((s) => [s.id, s.provenance])).toEqual([
      ['semtype:BinaryImage', 'ast'], ['semtype:DemoId', 'registry'], ['semtype:Macromolecule', 'annotation'], ['semtype:Molecule', 'annotation'], ['semtype:Text', 'ast'],
      ['semtype:demo-country', 'ast'], ['semtype:flag', 'ast'],
    ]);
  });

  it('resolves by-name calls case-insensitively and space-insensitively per source file with a count; unresolved targets are counted', async () => {
    const {rows, problems, manifest} = await graph;
    expect(rows('edges/calls').map((e) => [e.from, e.to, e.count])).toEqual([
      ['file:public/packages/Demo/detectors.js', 'func:Demo:detectMolecules', 2],
      ['file:public/packages/Demo/src/utils.ts', 'func:Demo:Calculate logD', 1],
      ['file:public/packages/Demo/src/utils.ts', 'func:Demo:Pareto Front', 1],
      ['file:public/packages/Demo/src/utils.ts', 'func:Demo:To HELM', 2],
      ['file:public/packages/Demo/src/utils.ts', 'func:Plain:helper', 1],
      ['file:public/packages/Plain/src/package.ts', 'func:Demo:toMolfile', 1],
    ]);
    expect(rows('edges/calls')[0]).toMatchObject({kind: 'by-name', derived_by: 'ast', confidence: 1, evidence: ['public/packages/Demo/detectors.js']});
    expect(problems.unresolved_ids).toContain('public/packages/Demo/src/utils.ts: call to Nowhere:missing names no registered function');
    expect(manifest.problems).toMatchObject({unresolved_ids: 3, dangling_edges: 0});
  });

  it('reads scripts with their language, environment reference, reference, sample and test; a file without a header is not a script', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/script'), 'func:Demo:Calculate logD')).toMatchObject({
      type: 'script', language: 'python', path: 'public/packages/Demo/scripts/sub/calc_logD.py', environment: 'env:Demo:demo-env', reference: 'https://en.wikipedia.org/wiki/Partition_coefficient',
      sample: 'chem/smiles.csv', test: true, tags: ['demo'], description: 'Calculates logD.', meta: {timeout: '900000', test: 'expect(1, 1) //cat: Types'},
      signature: '(dataframe table, column molecules: Molecule, double pH) -> dataframe', provenance: 'annotation',
    });
    expect(byId(rows('nodes/script'), 'func:Demo:Inline Env')).toMatchObject({language: 'python', test: false, meta: {environment: 'channels: [conda-forge], dependencies: [python=3.12, {pip: [cobra]}]'}});
    expect(byId(rows('nodes/script'), 'func:Demo:RunJs')).toMatchObject({language: 'js'});
    expect(rows('nodes/script')).toHaveLength(3);
    expect(rows('edges/environment')).toEqual([expect.objectContaining({type: 'ref', name: 'environment', from: 'func:Demo:Calculate logD', to: 'env:Demo:demo-env'})]);
    expect(byId(rows('nodes/script-environment'), 'env:Demo:demo-env')).toMatchObject({language: 'python', packages: ['python', 'pip', 'rdkit', 'pichemist'], path: 'public/packages/Demo/environments/demo-env.yaml', provenance: 'registry'});
  });

  it('reads every query of a .sql file with its connection reference; one without --connection: is skipped and makes the source partial', async () => {
    const {rows, problems, manifest} = await graph;
    expect(byId(rows('nodes/query'), 'func:Demo:compounds for @target')).toMatchObject({
      type: 'query', language: 'sql', line: 1, connection: 'conn:Demo:Demo', friendly_name: 'Browse | Compounds For Target', test: true, expected_rows: 5, cache: 'server',
      input_types: ['string'], meta: {test: "Dbtests:expectTable(CompoundsForTarget(), OpenFile('x.d42'))"},
    });
    expect(byId(rows('nodes/query'), 'func:Demo:users')).toMatchObject({line: 12, connection: 'conn:System:Datagrok', test: false});
    expect(rows('edges/connection').map((e) => [e.from, e.to])).toEqual([['func:Demo:compounds for @target', 'conn:Demo:Demo'], ['func:Demo:users', 'conn:System:Datagrok'],
      ['query:Demo:Demo App', 'conn:Demo:Demo']]);
    expect(byId(rows('nodes/connection'), 'conn:System:Datagrok')).toMatchObject({status: 'proposed'});
    expect(problems.invalid_rows).toEqual(["public/packages/Demo/queries/q.sql:19: query 'no connection' has no --connection:"]);
    expect(manifest.sources['ts-functions']).toBe('partial(2 rejected)');
  });

  it('reads connections without credentials', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/connection'), 'conn:Demo:Demo')).toEqual({
      id: 'conn:Demo:Demo', type: 'connection', name: 'Demo', batch: expect.any(String), db: 'demo', description: 'Demo DB', language: 'other', package: 'pkg:Demo',
      path: 'public/packages/Demo/connections/demo.json', provenance: 'registry', provider: 'Postgres', server: '${Demo<DockerContainer>}', source_layer: 'public', status: 'active', visibility: 'public',
    });
  });
});

describe('ts-functions extractor, the defects of the WO-3a review', () => {
  it('D1 reads package.g.ts, package.ts and package-test.ts: the generated file takes a name from package.ts silently, any other pair is shadowed', async () => {
    const {rows, problems} = await graph;
    expect(byId(rows('nodes/function'), 'func:Demo:info')).toMatchObject({path: 'public/packages/Demo/src/package.ts', line: 7});
    expect(byId(rows('nodes/app'), 'func:Demo:Demo App')).toMatchObject({path: 'public/packages/Demo/src/package.g.ts'});
    expect(problems.shadowed_headers).toEqual(["public/packages/Plain/src/package-test.ts:11: 'helper' is already registered in public/packages/Plain/src/package.ts"]);
  });

  it('D2 keeps the first of two declarations under one id, counts the collision apart from a merge, and gives a colliding query its own id', async () => {
    const {rows, problems} = await graph;
    expect(byId(rows('nodes/panel'), 'func:Demo:Dual')).toMatchObject({line: 153});
    expect(byId(rows('nodes/function'), 'func:Demo:Twice')).toMatchObject({line: 139, input_types: ['int']});
    expect(byId(rows('nodes/query'), 'query:Demo:Demo App')).toMatchObject({language: 'sql', connection: 'conn:Demo:Demo'});
    expect(byId(rows('nodes/app'), 'func:Demo:Demo App')).toMatchObject({language: 'ts'});
    // two registrations competing for one name are a collision, not the cross-extractor merge duplicate_ids counts
    expect(problems.registration_collisions).toEqual([
      'func:Demo:Twice: function at public/packages/Demo/src/package.g.ts:145 ignored; function at public/packages/Demo/src/package.g.ts:139 kept',
      'func:Demo:Dual: app at public/packages/Demo/src/package.g.ts:162 ignored; panel at public/packages/Demo/src/package.g.ts:153 kept',
      "public/packages/Demo/queries/q.sql: query 'Demo App' collides with the function of the same name in public/packages/Demo/src/package.g.ts; kept as query:Demo:Demo App",
    ]);
    expect(problems.duplicate_ids).toBeUndefined();
  });

  it('D2 re-checks a merged row against the winning type and demotes it to invalid', () => {
    const emitter = new Emitter(loadTypeSystem(path.join(fixture, KG_DIR)), 'b-test');
    const at = {name: 'Merged', path: 'public/packages/Demo/queries/q.sql', line: 1, package: 'pkg:Demo', provenance: 'annotation', source_layer: 'public'};
    emitter.node({type: 'function', id: 'func:Demo:Merged', language: 'ts', ...at});
    emitter.node({type: 'query', id: 'func:Demo:Merged', language: 'sql', connection: 'conn:Demo:Demo', ...at});
    const graph = emitter.finalize();
    expect(graph.nodes.map((n) => n.id)).not.toContain('func:Demo:Merged');
    expect(graph.invalid).toEqual([expect.objectContaining({id: 'func:Demo:Merged', type: 'query', problems: [expect.stringContaining('language')]})]);
  });

  it('D3 does not let a blank line inside a header end it, in all three parsers', async () => {
    const {rows} = await graph;
    expect(parseHeaderLines(['//name: gap', '', '//top-menu: Demo | Gap'], 3)).toMatchObject({line: 3, name: 'gap', keys: {'top-menu': ['Demo | Gap']}});
    expect(byId(rows('nodes/function'), 'func:Demo:Spaced Header')).toMatchObject({top_menu: 'Demo | Spaced'});
    expect(byId(rows('nodes/query'), 'func:Demo:users')).toMatchObject({line: 12, connection: 'conn:System:Datagrok'});
    expect(byId(rows('nodes/script'), 'func:Demo:Inline Env')).toMatchObject({language: 'python'});
  });

  it('D4 reads the three container layouts: a folder per Dockerfile, a single dockerfiles/Dockerfile, a published image in container.json', async () => {
    const {rows} = await graph;
    expect(rows('nodes/container').map((c) => [c.id, c.name, c.path])).toEqual([
      ['container:Demo:Demo', 'Demo', 'public/packages/Demo/dockerfiles/Dockerfile'],
      ['container:Demo:demo', 'demo', 'public/packages/Demo/dockerfiles/demo/Dockerfile'],
      ['container:Demo:published', 'published', 'public/packages/Demo/dockerfiles/published/container.json'],
    ]);
    expect(rows('nodes/container')[0].base).toBeUndefined();
  });

  it('D5 resolves a call by the exact name before the spaceless one, and reports a collision instead of guessing', async () => {
    const {rows, problems} = await graph;
    expect(rows('edges/calls').find((e) => e.to.startsWith('func:Demo:Pareto'))).toMatchObject({to: 'func:Demo:Pareto Front'});
    expect(rows('edges/calls').some((e) => e.to === 'func:Demo:Trade Off' || e.to === 'func:Demo:Trade off')).toBe(false);
    expect(problems.ambiguous_calls).toEqual(['public/packages/Demo/src/utils.ts: call to Demo:TradeOff matches func:Demo:Trade Off, func:Demo:Trade off']);
  });

  it('D6 reads a JavaScript entry point at the package root and under src', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/function'), 'func:Demo:rootJs')).toMatchObject({path: 'public/packages/Demo/package.js', language: 'js', description: 'A JavaScript entry point at the package root'});
    expect(byId(rows('nodes/function'), 'func:Plain:srcJs')).toMatchObject({path: 'public/packages/Plain/src/package.js', language: 'js'});
  });

  it('D7 reads columnless from meta.columnlessFilter and leaves an editor without edits, so //editor-for: stays in meta', async () => {
    const {rows} = await graph;
    expect(byId(rows('nodes/filter'), 'func:Demo:Substructure Filter')).toMatchObject({columnless: true, primary: true});
    expect(byId(rows('nodes/editor'), 'func:Demo:Column Editor')).toMatchObject({meta: {'editor-for': 'AddNewColumn'}});
    expect(byId(rows('nodes/editor'), 'func:Demo:Column Editor').edits).toBeUndefined();
  });

  it('D8 detects the semantic type of a ternary return and of a constant map, and counts a constant it cannot resolve', async () => {
    const {rows, problems} = await graph;
    const detects = rows('edges/targets-semtype').filter((e) => e.role === 'detects').map((e) => [e.from, e.to]);
    expect(detects).toContainEqual(['func:Demo:detectFlags', 'semtype:flag']);
    expect(detects).toContainEqual(['func:Demo:detectCountries', 'semtype:demo-country']);
    expect(detects.some((d) => d[0] === 'func:Demo:detectNowhere')).toBe(false);
    expect(problems.unresolved_ids).toContain('public/packages/Demo/detectors.js: ELSEWHERE_SEMTYPES.NOPE names no semantic type this file declares');
  });
});
