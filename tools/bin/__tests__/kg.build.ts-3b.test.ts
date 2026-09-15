/// `grok kg build` WO-3b (build-plan.md): the ts-declarations, ts-imports and ts-uses extractors against the
/// mini monorepo under fixtures/kg/build, with ts-packages and ts-functions for the owners they attach to.
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {apiTokens} from '../utils/kg/build/extract/ts/uses';
import {currentDir} from '../utils/kg/build/write';
import {kg} from '../commands/kg';

const fixture = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg', 'build');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const API = 'public/js-api';
const VIEWER = 'public/packages/Demo/src/viewer.ts';

async function build(): Promise<{manifest: any, rows: (file: string) => any[], problems: Record<string, string[]>}> {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-3b-'));
  fs.cpSync(fixture, repo, {recursive: true});
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    await kg({_: ['kg', 'build'], kg: path.join(repo, KG_DIR), only: 'ts-packages,ts-functions,ts-declarations,ts-imports,ts-uses', db: false, output: 'json'});
    expect(error.mock.calls).toEqual([]);
    expect(process.exitCode).toBeUndefined();
    const out = currentDir(path.join(repo, '.kg'))!;
    const rows = (file: string) => {
      const p = path.join(out, `data/${file}.jsonl`);
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
const pairs = (rows: any[], from: string) => rows.filter((e) => e.from === from).map((e) => e.to);

describe('ts-declarations extractor (build-plan.md WO-3b)', () => {
  it('walks js-api, packages and libraries into source-file nodes with loc, language and generated, declared by their owner', async () => {
    const {rows, manifest} = await graph;
    expect(manifest.sources).toMatchObject({'ts-declarations': 'ok', 'ts-imports': 'ok', 'ts-uses': 'ok'});
    expect(byId(rows('nodes/source-file'), `file:${VIEWER}`)).toMatchObject({name: 'viewer.ts', path: VIEWER, loc: 38, language: 'ts', generated: false, package: 'pkg:Demo', provenance: 'filesystem', source_layer: 'public', visibility: 'public'});
    expect(byId(rows('nodes/source-file'), `file:${API}/src/api/ddt.api.g.ts`)).toMatchObject({generated: true});
    expect(byId(rows('nodes/source-file'), `file:${API}/ui.ts`).package).toBeUndefined();
    const declares = rows('edges/declares');
    expect(declares.find((e) => e.from === 'lib:js-api' && e.to === `file:${API}/ui.ts`)).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [`${API}/ui.ts`]});
    expect(pairs(declares, 'lib:utils')).toContain('file:public/libraries/utils/src/test.ts');
    expect(pairs(declares, 'pkg:Demo')).toContain(`file:${VIEWER}`);
    expect(rows('nodes/source-file').some((f) => f.path.endsWith('.d.ts') || f.path.endsWith('.css'))).toBe(false);
  });

  it('emits declarations with kind, exported, documented, signature and line; ids per ids.ts', async () => {
    const {rows} = await graph;
    const decls = rows('nodes/declaration');
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame`)).toEqual({
      id: `decl:${API}/src/dataframe.ts#DataFrame`, type: 'declaration', name: 'DataFrame', batch: expect.any(String), deprecated: false, documented: true, exported: true,
      generated: false, kind: 'class', language: 'ts', line: 2, path: `${API}/src/dataframe.ts`, provenance: 'ast', public_api: true, signature: 'export class DataFrame',
      source_layer: 'public', status: 'active', visibility: 'public',
    });
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.fromCsv`)).toMatchObject({name: 'fromCsv', kind: 'method', exported: true, public_api: true,
      documented: false, line: 23, signature: 'static fromCsv(csv: string, options?: {delimiter?: string}): DataFrame'});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.fromCsv`)).not.toHaveProperty('container');
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.hidden`)).toMatchObject({kind: 'method', exported: false, public_api: false});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame._name`)).toMatchObject({kind: 'prop', exported: false, signature: "private _name: string = ''"});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#IDisposable.dispose`)).toMatchObject({kind: 'method'});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#Internal`)).toMatchObject({kind: 'class', exported: false, public_api: false});
    expect(decls.filter((d) => d.path === `${API}/src/dataframe.ts`).map((d) => [d.name, d.kind]).sort()).toEqual([
      ['DataFrame', 'class'], ['IDisposable', 'interface'], ['Internal', 'class'], ['LogLevel', 'enum'], ['Predicate', 'type'], ['_name', 'prop'], ['dispose', 'method'],
      ['fromCsv', 'method'], ['fromText', 'method'], ['hidden', 'method'], ['m', 'method'], ['name', 'getter'], ['name', 'setter'], ['rowCount', 'getter'],
    ]);
    expect(byId(decls, `decl:${API}/src/const.ts#SEMTYPE`)).toMatchObject({kind: 'const', documented: true, signature: 'export const SEMTYPE'});
    expect(byId(decls, 'decl:public/packages/Plain/src/package.ts#plainApp')).toMatchObject({kind: 'const', signature: 'export const plainApp = async (): Promise<void> =>'});
  });

  it('gives a getter and a setter of one name the :get and :set ids, and overloads one id with the merged signature', async () => {
    const {rows} = await graph;
    const decls = rows('nodes/declaration');
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.name:get`)).toMatchObject({name: 'name', kind: 'getter', documented: true, line: 6, signature: 'get name(): string'});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.name:set`)).toMatchObject({name: 'name', kind: 'setter', documented: false, line: 10, signature: 'set name(s: string)'});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.name`)).toBeUndefined();
    const button = decls.filter((d) => d.path === `${API}/ui.ts` && d.name === 'button');
    expect(button).toHaveLength(1);
    expect(button[0]).toMatchObject({id: `decl:${API}/ui.ts#button`, kind: 'function', line: 8});
    expect(button[0].signature).toMatch(/^export function button\(text: string\): HTMLButtonElement \| export function button\(text: string, onClick: \(\) => void\): HTMLButtonElement/);
    expect(button[0].signature.length).toBe(160);
  });

  it('caps the signature at 160 characters', async () => {
    const {rows} = await graph;
    const long = byId(rows('nodes/declaration'), `decl:${API}/ui.ts#longSignature`);
    expect(long.signature).toHaveLength(160);
    expect(long.signature.startsWith('export function longSignature(aVeryLongParameterName: string, anotherVeryLongParameterName: number,')).toBe(true);
  });

  it('sets public_api only for exported js-api declarations, and deprecated from a JSDoc @deprecated tag', async () => {
    const {rows} = await graph;
    const decls = rows('nodes/declaration');
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.fromText`)).toMatchObject({deprecated: true, documented: true, public_api: true});
    expect(byId(decls, `decl:${API}/src/dataframe.ts#DataFrame.fromCsv`)).toMatchObject({deprecated: false});
    expect(byId(decls, `decl:${API}/ui.ts#internal`)).toMatchObject({exported: false, public_api: false});
    expect(byId(decls, 'decl:public/libraries/utils/src/test.ts#something')).toMatchObject({exported: true, public_api: false});
    expect(byId(decls, 'decl:public/packages/Demo/src/base.ts#Base')).toMatchObject({exported: true, public_api: false});
    expect(decls.filter((d) => d.public_api).every((d) => d.path.startsWith(`${API}/`) && d.exported)).toBe(true);
  });

  it('marks declarations of *.g.ts generated and keeps only their public_api members', async () => {
    const {rows} = await graph;
    const generated = rows('nodes/declaration').filter((d) => d.path === `${API}/src/api/ddt.api.g.ts`);
    expect(generated.map((d) => [d.name, d.kind, d.generated, d.public_api]).sort()).toEqual([
      ['DataFrameApi', 'class', true, true], ['Hidden', 'class', true, false], ['generated', 'method', true, true],
    ]);
  });

  it('declares each declaration from its file, each member from its class; a namespace is a const holding its members', async () => {
    const {rows} = await graph;
    const declares = rows('edges/declares');
    const df = `decl:${API}/src/dataframe.ts#DataFrame`;
    expect(pairs(declares, `file:${API}/src/dataframe.ts`).sort()).toEqual([df, `decl:${API}/src/dataframe.ts#IDisposable`, `decl:${API}/src/dataframe.ts#Internal`, `decl:${API}/src/dataframe.ts#LogLevel`, `decl:${API}/src/dataframe.ts#Predicate`]);
    expect(pairs(declares, df).sort()).toEqual([`${df}._name`, `${df}.fromCsv`, `${df}.fromText`, `${df}.hidden`, `${df}.name:get`, `${df}.name:set`, `${df}.rowCount:get`]);
    expect(declares.find((e) => e.from === df && e.to === `${df}.fromCsv`)).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [`${API}/src/dataframe.ts`]});
    expect(declares.some((e) => e.from === `file:${API}/src/dataframe.ts` && e.to === `${df}.fromCsv`)).toBe(false);
    expect(byId(rows('nodes/declaration'), `decl:${API}/ui.ts#input`)).toMatchObject({kind: 'const', signature: 'export namespace input'});
    expect(byId(rows('nodes/declaration'), `decl:${API}/ui.ts#input.string`)).toMatchObject({kind: 'function', public_api: true});
    expect(pairs(declares, `decl:${API}/ui.ts#input`)).toEqual([`decl:${API}/ui.ts#input.string`]);
    expect(rows('edges/container')).toEqual([]);
  });

  it('resolves extends and implements in the same file, through an import, and across files of the same package; a stranger is dropped and counted', async () => {
    const {rows, problems} = await graph;
    expect(rows('edges/extends').map((e) => [e.from, e.to])).toEqual([
      [`decl:${API}/src/chem.ts#chem.SearchOptions`, `decl:${API}/src/chem.ts#chem.Options`],
      [`decl:${API}/src/viewer.ts#JsViewer`, `decl:${API}/src/viewer.ts#Viewer`],
      [`decl:${VIEWER}#Child`, 'decl:public/packages/Demo/src/base.ts#Base'],
      [`decl:${VIEWER}#DemoViewer`, `decl:${API}/src/viewer.ts#JsViewer`],
    ]);
    expect(rows('edges/extends')[3]).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [VIEWER]});
    expect(rows('edges/implements').map((e) => [e.from, e.to])).toEqual([[`decl:${VIEWER}#DemoViewer`, `decl:${API}/src/dataframe.ts#IDisposable`]]);
    expect(problems.unresolved_ids).toContain(`${VIEWER}:37: Lost extends Missing, which no file in reach declares`);
    expect(byId(rows('nodes/declaration'), `decl:${VIEWER}#Lost`)).toMatchObject({kind: 'class', signature: 'export class Lost extends Missing'});
  });
});

describe('ts-imports extractor (build-plan.md WO-3b)', () => {
  it('resolves every specifier form to a file, a library or a package with the symbols named; bare names and assets stay out', async () => {
    const {rows, problems} = await graph;
    const imports = rows('edges/imports').filter((e) => e.from === `file:${VIEWER}`).map((e) => [e.to, e.symbols]);
    expect(imports).toEqual([
      ['file:public/libraries/utils/src/test.ts', ['other', 'something']],
      ['file:public/packages/Demo/src/base.ts', ['Base']],
      ['file:public/packages/Demo/src/common/index.ts', ['shared']],
      ['file:public/packages/Demo/src/utils.ts', ['callThings']],
      ['file:public/packages/Plain/src/package.ts', ['helper']],
      ['lib:js-api', ['*']],
      ['lib:utils', undefined],
    ]);
    expect(rows('edges/imports').find((e) => e.from === `file:${VIEWER}` && e.to === 'lib:js-api')).toMatchObject({derived_by: 'ast', confidence: 1, evidence: [VIEWER]});
    expect(problems.unresolved_ids).toContain(`${VIEWER}: import './missing' resolves to no file, library or package`);
    expect(problems.unresolved_ids.some((p) => p.includes("'rxjs'") || p.includes("'./styles.css'"))).toBe(false);
  });

  it('follows re-exports and named imports inside js-api to their files', async () => {
    const {rows} = await graph;
    expect(rows('edges/imports').filter((e) => e.from === `file:${API}/grok.ts`).map((e) => [e.to, e.symbols])).toEqual([
      [`file:${API}/src/chem.ts`, ['*']], [`file:${API}/src/logger.ts`, ['Logger']], [`file:${API}/src/shell.ts`, ['Shell']],
    ]);
    expect(pairs(rows('edges/imports'), `file:${API}/dg.ts`)).toEqual([`file:${API}/src/const.ts`, `file:${API}/src/dataframe.ts`, `file:${API}/src/shell.ts`, `file:${API}/src/u2core/index.ts`, `file:${API}/src/viewer.ts`]);
  });
});

describe('ts-uses extractor (build-plan.md WO-3b)', () => {
  it('counts qualified DG, ui and grok tokens outside comments', () => {
    expect([...apiTokens("// DG.Column\n/* ui.span */ DG.DataFrame.fromCsv('x'); grok.shell.info(DG.DataFrame); ui.div(); const s = 'grok.shell.error'; grokery.x; DG.a.b.c.d")]).toEqual([
      ['DG.DataFrame.fromCsv', 1], ['grok.shell.info', 1], ['DG.DataFrame', 1], ['ui.div', 1], ['grok.shell.error', 1], ['DG.a.b.c', 1],
    ]);
  });

  it('resolves DG.X to the class, enum, const or type (DG.U2.X through the namespace re-export), ui.x to ui.ts, grok.a.b to the member of the class behind a, with kind and count', async () => {
    const {rows} = await graph;
    const uses = rows('edges/uses').filter((e) => e.from === `file:${VIEWER}`).map((e) => [e.to, e.kind, e.count]);
    expect(uses).toEqual([
      [`decl:${API}/src/chem.ts#chem.similarity`, 'function', 1],
      [`decl:${API}/src/const.ts#SEMTYPE`, 'enum', 1],
      [`decl:${API}/src/dataframe.ts#DataFrame`, 'class', 1],
      [`decl:${API}/src/dataframe.ts#IDisposable`, 'type', 1],
      [`decl:${API}/src/dataframe.ts#LogLevel`, 'enum', 1],
      [`decl:${API}/src/logger.ts#Logger.info`, 'function', 1],
      [`decl:${API}/src/shell.ts#Shell.info`, 'function', 2],
      [`decl:${API}/src/u2core/index.ts#Control`, 'class', 1],
      [`decl:${API}/src/viewer.ts#JsViewer`, 'class', 1],
      [`decl:${API}/ui.ts#div`, 'ui', 1],
      [`decl:${API}/ui.ts#input.string`, 'ui', 1],
    ]);
    // a token scan cannot tell a binding from a shadowed name or a string, so the edge claims less than certainty
    expect(rows('edges/uses').find((e) => e.from === `file:${VIEWER}` && e.to === `decl:${API}/src/shell.ts#Shell.info`)).toMatchObject({derived_by: 'ast', confidence: 0.8, evidence: [VIEWER]});
    expect(rows('edges/uses').filter((e) => e.from === 'file:public/packages/Demo/src/package.g.ts').map((e) => e.to)).toEqual([`decl:${API}/src/dataframe.ts#DataFrame`]);
  });

  it('counts unresolved tokens once per token in the manifest problems, with the number of files', async () => {
    const {problems, manifest} = await graph;
    expect(problems.unresolved_ids).toContain(`uses: DG.Nowhere.x names no exported JS API declaration (1 file, first ${VIEWER})`);
    expect(problems.unresolved_ids).toContain(`uses: grok.nowhere.y names no exported JS API declaration (1 file, first ${VIEWER})`);
    expect(problems.unresolved_ids).toContain('uses: grok.functions.call names no exported JS API declaration (2 files, first public/packages/Demo/src/utils.ts)');
    expect(problems.unresolved_ids.filter((p) => p.startsWith('uses: DG.Nowhere.x'))).toHaveLength(1);
    expect(manifest.problems.unresolved_ids).toBe(problems.unresolved_ids.length);
    expect(manifest.problems.dangling_edges).toBe(0);
  });
});
