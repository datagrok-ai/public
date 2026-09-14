/// Registered functions of every package (build-plan.md WO-3a): annotation headers in every entry file
/// (`src/package.g.ts`, `src/package.ts`, `src/package-test.ts`, `package.js`) and `detectors.js`, scripts,
/// queries, connections, environments and containers as their own node types, `targets-semtype` per role,
/// and by-name `calls` per source file.
import * as fs from 'fs';
import * as path from 'path';
import * as yaml from 'js-yaml';
import {globSync} from 'glob';
import {FUNC_TYPES} from '../../../../const';
import {Emitter} from '../../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../../registry';
import {pkgId, fileId, declId, funcId, connId, envId, containerId, semtypeId, languageOf} from '../../ids';
import {Header, HeaderBlock, parseFunctionHeaders, parseScriptHeader, parseQueryHeaders} from '../../annotations';
import {countLines} from '../homes';
import {PackageFolder, listPackages} from './packages';

/** function.yaml: the subtype of a function with several roles is that of its highest-precedence role. */
const ROLE_PRECEDENCE = ['app', 'viewer', 'filter', 'cellRenderer', 'semTypeDetector', 'fileHandler', 'fileExporter', 'fileViewer', 'scriptHandler', 'panel', 'editor', 'valueEditor', 'cellEditor', 'init', 'autostart'];
const ROLE_TYPES: Record<string, string> = {
  app: 'app', viewer: 'viewer', filter: 'filter', cellRenderer: 'cell-renderer', semTypeDetector: 'sem-type-detector', fileHandler: 'file-handler',
  fileExporter: 'file-handler', fileViewer: 'file-viewer', scriptHandler: 'script-handler', panel: 'panel', editor: 'editor', valueEditor: 'editor',
  cellEditor: 'editor', init: 'lifecycle-hook', autostart: 'lifecycle-hook',
};
const KNOWN_ROLES = new Set<string>(Object.values(FUNC_TYPES));
const CACHE: Record<string, string> = {all: 'all', server: 'server', client: 'client', true: 'all'};
const SCRIPT_LANGUAGES: Record<string, string> = {javascript: 'js', nodejs: 'js', pyodide: 'python', python: 'python', r: 'r', julia: 'julia', octave: 'octave', grok: 'grok'};
const SCRIPT_LANGUAGE_ENUM = ['python', 'r', 'julia', 'octave', 'js', 'grok'];
const SCRIPT_GLOB = 'scripts/**/*.{py,R,r,jl,m,js,grok}';
const SOURCE_IGNORE = ['**/node_modules/**', '**/dist/**'];
/** Its scripts/ folder holds API samples (WO-3c), not scripts. */
const SAMPLES_PACKAGE = 'ApiSamples';
const URL = /^[a-z][a-z0-9+.-]*:\/\//i;
const CALL_PATTERNS = [
  /grok\.functions\.(?:call|eval)\(\s*['"`](\w+):(\w+)/g,
  /DG\.Func\.byName\(\s*['"`](\w+):(\w+)/g,
  /DG\.Func\.find\(\s*\{\s*package\s*:\s*['"`](\w+)['"`]\s*,\s*name\s*:\s*['"`](\w+)/g,
];
/** Header keys that become members or edges; every other key lands in `meta` verbatim. */
const LIFTED_KEYS = new Set(['name', 'description', 'input', 'output', 'tags', 'friendlyName', 'top-menu', 'feature']);
/** Where a package registers functions with headers, in precedence order: the generated file is the decorator-lowered
 * form of what `src/package.ts` declares, so it wins over it in silence. */
const ENTRY_FILES = ['src/package.g.ts', 'src/package.ts', 'src/package.js', 'package.js', 'src/package-test.ts', 'detectors.js'];
const GENERATED_ENTRY = 'src/package.g.ts';
const LOWERED_ENTRY = 'src/package.ts';
/** What a detector returns or assigns, plainly and as the `cond ? 'x' : null` a detector may end with. */
const DETECTED = [
  /(?:return|\.semType\s*=)\s+(?:DG\.SEMTYPE\.([A-Z][A-Z0-9_]*)|(['"])([^'"\n]+)\2|([A-Z][A-Z0-9_]{2,})(?:\.([A-Za-z_][\w]*))?)/g,
  /(?:return|\.semType\s*=)[^;\n]*\?\s*(?:DG\.SEMTYPE\.([A-Z][A-Z0-9_]*)|(['"])([^'"\n]+)\2|([A-Z][A-Z0-9_]{2,})(?:\.([A-Za-z_][\w]*))?)\s*:\s*null/g,
];

interface SemtypeUse {
  role: string;
  semtype: string;
  derived_by: string;
  confidence: number;
}

interface FunctionSpec {
  id: string;
  name: string;
  language: string;
  path: string;
  line?: number;
  pkg: PackageFolder;
  /** `function` lets the roles pick a subtype; a script or a query keeps its type. */
  type: 'function' | 'script' | 'query';
}

export const functionsExtractor: Extractor = {
  name: 'ts-functions',
  layer: 'public',
  modes: ['full'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const layer = new FunctionLayer(ctx.repoRoot, emitter);
    const packages = listPackages(ctx.repoRoot);
    for (const pkg of packages) layer.emitPackage(pkg);
    for (const pkg of packages) layer.emitCalls(pkg);
    emitter.source('ts-functions', layer.skipped ? 'partial' : 'ok');
  },
};

class FunctionLayer {
  skipped = 0;
  /** `<pkg>:<name>` lowercased -> the functions registered under it; [exact] keeps the spaces of the name, [loose] strips them. */
  private exact = new Map<string, Set<string>>();
  private loose = new Map<string, Set<string>>();
  /** The names the entry files of the package at hand registered, and where. */
  private declared = new Map<string, string>();
  private files = new Set<string>();
  private semtypes: Record<string, string>;

  constructor(private repoRoot: string, private emitter: Emitter) {
    this.semtypes = loadSemtypes(repoRoot);
  }

  emitPackage(pkg: PackageFolder): void {
    this.declared.clear();
    for (const entry of ENTRY_FILES) {
      const file = `${pkg.dir}/${entry}`;
      if (this.exists(file)) this.emitFunctionFile(pkg, file, entry);
    }
    if (pkg.folder !== SAMPLES_PACKAGE)
      for (const file of this.glob(pkg, SCRIPT_GLOB)) this.emitScript(pkg, file);
    for (const file of this.glob(pkg, 'queries/**/*.sql')) this.emitQueries(pkg, file);
    for (const file of this.glob(pkg, 'connections/*.json')) this.emitConnection(pkg, file);
    for (const file of this.glob(pkg, 'environments/*.{yaml,yml}')) this.emitEnvironment(pkg, file);
    this.emitContainers(pkg);
  }

  /** `calls` by name per source file, resolved case-insensitively against every registered name; unresolved and ambiguous targets are counted. */
  emitCalls(pkg: PackageFolder): void {
    const files = [...this.glob(pkg, 'src/**/*.{ts,js}').filter((f) => !f.endsWith('.d.ts')), `${pkg.dir}/detectors.js`].filter((f) => this.exists(f));
    for (const file of files) {
      const text = this.read(file);
      const counts = new Map<string, number>();
      for (const pattern of CALL_PATTERNS)
        for (const m of text.matchAll(pattern)) {
          const candidates = this.exact.get(exactKey(m[1], m[2])) ?? this.loose.get(looseKey(m[1], m[2]));
          if (!candidates)
            this.emitter.problem('unresolved_ids', `${file}: call to ${m[1]}:${m[2]} names no registered function`);
          else if (candidates.size > 1)
            this.emitter.problem('ambiguous_calls', `${file}: call to ${m[1]}:${m[2]} matches ${[...candidates].sort().join(', ')}`);
          else {
            const to = [...candidates][0];
            counts.set(to, (counts.get(to) ?? 0) + 1);
          }
        }
      if (!counts.size) continue;
      this.fileNode(pkg, file);
      for (const [to, count] of counts)
        this.emitter.edge({type: 'calls', from: fileId(file), to, kind: 'by-name', count, derived_by: 'ast', confidence: 1, evidence: [file]});
    }
  }

  private emitFunctionFile(pkg: PackageFolder, file: string, entry: string): void {
    const text = this.read(file);
    this.fileNode(pkg, file);
    for (const block of parseFunctionHeaders(text)) {
      const name = block.header.name ?? block.declaration.name;
      const first = this.declared.get(name);
      if (first !== undefined && first !== entry) {
        if (!(first === GENERATED_ENTRY && entry === LOWERED_ENTRY))
          this.emitter.problem('shadowed_headers', `${file}:${block.declaration.line}: '${name}' is already registered in ${pkg.dir}/${first}`);
        continue;
      }
      this.declared.set(name, entry);
      const id = funcId(pkg.folder, name);
      const {row, uses} = this.functionRow(block.header, {id, name, language: languageOf(file), path: file, line: block.declaration.line, pkg, type: 'function'});
      if (row.type === 'sem-type-detector' && !uses.some((u) => u.role === 'detects'))
        for (const semtype of this.detectedSemtypes(block, text, file)) uses.push({role: 'detects', semtype, derived_by: 'ast', confidence: 0.9});
      this.emitFunction(row, uses, block.header, pkg, file);
      this.emitter.edge({type: 'declares', from: fileId(file), to: id, derived_by: 'ast', confidence: 1, evidence: [file]});
      if (block.declaration.kind === 'method' && block.declaration.owner) {
        const owner = declId(file, block.declaration.owner);
        this.emitter.stub(owner, 'declaration', block.declaration.owner, 'ast', {language: languageOf(file), path: file, kind: 'class'});
        this.emitter.edge({type: 'declares', from: owner, to: id, derived_by: 'ast', confidence: 1, evidence: [file]});
      }
      this.register(pkg.folder, name, id);
      this.register(pkg.folder, block.declaration.name, id);
    }
  }

  private emitScript(pkg: PackageFolder, file: string): void {
    const text = this.read(file);
    const header = parseScriptHeader(text, languageOf(file));
    if (!header || (!header.name && !header.keys.language)) return;
    const tag = header.keys.language?.[0]?.toLowerCase();
    const byTag = tag === undefined ? undefined : SCRIPT_LANGUAGES[tag];
    const byExt = languageOf(file);
    const language = byTag ?? (SCRIPT_LANGUAGE_ENUM.includes(byExt) ? byExt : 'other');
    const name = header.name ?? path.posix.basename(file).replace(/\.[^.]+$/, '');
    const {row, uses} = this.functionRow(header, {id: this.ownId(pkg, 'script', name, file), name, language, path: file, pkg, type: 'script'});
    const meta = row.meta as Record<string, string>;
    delete meta.language;
    const environment = header.keys.environment?.[0];
    if (environment && /^[\w.-]+$/.test(environment)) {
      row.environment = envId(pkg.folder, environment);
      delete meta.environment;
    }
    const reference = header.keys.reference?.[0];
    if (reference && URL.test(reference)) {
      row.reference = reference;
      delete meta.reference;
    }
    if (header.keys.sample?.[0]) {
      row.sample = header.keys.sample[0];
      delete meta.sample;
    }
    if (header.keys.test) row.test = true;
    this.emitFunction(row, uses, header, pkg, file);
    this.register(pkg.folder, name, row.id as string);
  }

  private emitQueries(pkg: PackageFolder, file: string): void {
    const text = this.read(file);
    for (const header of parseQueryHeaders(text)) {
      const name = header.name!;
      const connection = header.keys.connection?.[0];
      if (!connection) {
        this.emitter.problem('invalid_rows', `${file}:${header.line}: query '${name}' has no --connection:`);
        this.skipped++;
        continue;
      }
      const {row, uses} = this.functionRow(header, {id: this.ownId(pkg, 'query', name, file), name, language: 'sql', path: file, line: header.line, pkg, type: 'query'});
      const meta = row.meta as Record<string, string>;
      const colon = connection.indexOf(':');
      row.connection = colon < 0 ? connId(pkg.folder, connection) : connId(connection.slice(0, colon), connection.slice(colon + 1));
      delete meta.connection;
      if (header.keys.test) row.test = true;
      const expected = Number(meta.testExpectedRows);
      if (meta.testExpectedRows !== undefined && Number.isFinite(expected)) {
        row.expected_rows = expected;
        delete meta.testExpectedRows;
      }
      this.emitFunction(row, uses, header, pkg, file);
      this.register(pkg.folder, name, row.id as string);
    }
  }

  /** `func:<Pkg>:<name>`, unless a header function of the package already answers to that name: a script and a query are
   * separate registrations and must not merge into it. */
  private ownId(pkg: PackageFolder, scheme: 'script' | 'query', name: string, file: string): string {
    if (!this.declared.has(name)) return funcId(pkg.folder, name);
    const id = `${scheme}:${pkg.folder}:${name}`;
    this.emitter.problem('duplicate_ids', `${file}: ${scheme} '${name}' collides with the function of the same name in ${pkg.dir}/${this.declared.get(name)}; kept as ${id}`);
    return id;
  }

  /** Provider and endpoint only: credentials never enter the graph. */
  private emitConnection(pkg: PackageFolder, file: string): void {
    let json: any;
    try {
      json = JSON.parse(this.read(file));
    } catch {
      this.emitter.problem('invalid_rows', `${file}: not valid JSON`);
      this.skipped++;
      return;
    }
    if (typeof json?.dataSource !== 'string') return;
    const name = typeof json.name === 'string' ? json.name : path.posix.basename(file, '.json');
    const id = connId(pkg.folder, name);
    this.emitter.node({type: 'connection', id, name, description: json.description, language: 'other', path: file, package: pkgId(pkg.folder), provider: json.dataSource,
      server: json.parameters?.server, db: json.parameters?.db, provenance: 'registry', source_layer: 'public'});
    this.declares(pkg, id, file);
  }

  private emitEnvironment(pkg: PackageFolder, file: string): void {
    let data: any;
    try {
      data = yaml.load(this.read(file));
    } catch {
      this.emitter.problem('invalid_rows', `${file}: not valid YAML`);
      this.skipped++;
      return;
    }
    const name = typeof data?.name === 'string' ? data.name : path.posix.basename(file).replace(/\.ya?ml$/, '');
    const id = envId(pkg.folder, name);
    const packages = dependencyNames(data?.dependencies);
    const language = packages.some((p) => /^python\b/.test(p)) ? 'python' : packages.some((p) => /^r(-|$)/.test(p)) ? 'r' : 'other';
    this.emitter.node({type: 'script-environment', id, name, language, path: file, package: pkgId(pkg.folder), packages: packages.length ? packages : undefined,
      provenance: 'registry', source_layer: 'public'});
    this.declares(pkg, id, file);
  }

  /** The three layouts `grok publish` builds from (publish.ts:55-106): a folder per container, one `dockerfiles/Dockerfile`
   * whose image is the package itself, and a folder naming an already published image in `container.json`. */
  private emitContainers(pkg: PackageFolder): void {
    for (const file of this.glob(pkg, 'dockerfiles/*/Dockerfile'))
      this.emitContainer(pkg, path.posix.basename(path.posix.dirname(file)), file);
    const single = `${pkg.dir}/dockerfiles/Dockerfile`;
    if (this.exists(single)) this.emitContainer(pkg, pkg.folder, single);
    for (const file of this.glob(pkg, 'dockerfiles/*/container.json'))
      if (!this.exists(`${path.posix.dirname(file)}/Dockerfile`)) this.emitContainer(pkg, path.posix.basename(path.posix.dirname(file)), file);
  }

  /** `base` is left out: it references an image node no extractor produces yet (build-plan.md WO-3a). */
  private emitContainer(pkg: PackageFolder, name: string, file: string): void {
    const id = containerId(pkg.folder, name);
    this.emitter.node({type: 'container', id, name, language: 'other', path: file, package: pkgId(pkg.folder), provenance: 'filesystem', source_layer: 'public'});
    this.declares(pkg, id, file);
  }

  /** The node row of a header: members lifted by name, the subtype by role precedence with its members, everything else in `meta`. */
  private functionRow(header: Header, spec: FunctionSpec): {row: Row, uses: SemtypeUse[]} {
    const meta = {...header.meta};
    const roles = [...new Set([...(meta.role ?? '').split(',').map((r) => r.trim()).filter(Boolean), ...header.tags.filter((t) => KNOWN_ROLES.has(t))])];
    delete meta.role;
    const row: Row = {type: spec.type, id: spec.id, name: spec.name, description: header.description, language: spec.language, path: spec.path, line: spec.line,
      package: pkgId(spec.pkg.folder), friendly_name: header.keys.friendlyName?.[0], top_menu: header.keys['top-menu']?.[0], roles: roles.length ? roles : undefined,
      tags: header.tags.length ? header.tags : undefined, provenance: 'annotation', source_layer: 'public'};
    if (header.inputs.length || header.outputs.length) {
      row.signature = `(${header.inputs.map((p) => `${p.type} ${p.name}${p.options.semType ? `: ${p.options.semType}` : ''}`).join(', ')}) -> ${header.outputs.map((p) => p.type).join(' | ') || 'void'}`;
      row.input_types = header.inputs.map((p) => p.type);
      row.output_types = header.outputs.map((p) => p.type);
    }
    const lifted = new Set(LIFTED_KEYS);
    const helpUrl = header.keys['help-url']?.[0];
    if (helpUrl && URL.test(helpUrl)) {
      row.help_url = helpUrl;
      lifted.add('help-url');
    }
    if (CACHE[meta.cache]) {
      row.cache = CACHE[meta.cache];
      delete meta.cache;
    }
    if (meta.demoPath !== undefined) {
      row.demo_path = meta.demoPath;
      delete meta.demoPath;
    }
    delete meta.feature;
    const uses: SemtypeUse[] = [];
    for (const p of header.inputs) if (p.options.semType) uses.push({role: 'consumes', semtype: p.options.semType, derived_by: 'annotation', confidence: 1});
    for (const p of header.outputs) if (p.options.semType) uses.push({role: 'produces', semtype: p.options.semType, derived_by: 'annotation', confidence: 1});
    if (spec.type === 'function') {
      const role = ROLE_PRECEDENCE.find((r) => roles.includes(r) || (r === 'fileHandler' && roles.includes('file-handler')));
      if (role) {
        const sub = this.subtypeMembers(ROLE_TYPES[role], roles, header, {...meta}, lifted);
        if (sub) {
          row.type = ROLE_TYPES[role];
          Object.assign(row, sub.members);
          for (const k of sub.consumed) delete meta[k];
          uses.push(...sub.uses);
        }
      }
    }
    for (const [key, values] of Object.entries(header.keys))
      if (!lifted.has(key) && !key.startsWith('meta.')) meta[key] = values.join('\n');
    row.meta = meta;
    return {row, uses};
  }

  /** The members of a role subtype from `meta`; null when a required member is missing, so the row stays a plain function. */
  private subtypeMembers(type: string, roles: string[], header: Header, meta: Record<string, string>, lifted: Set<string>):
    {members: Row, consumed: string[], uses: SemtypeUse[]} | null {
    const consumed: string[] = [];
    const uses: SemtypeUse[] = [];
    const lift = (key: string) => {
      consumed.push(key);
      return meta[key];
    };
    const flag = (key: string) => meta[key] === undefined ? undefined : lift(key) === 'true';
    const list = (key: string) => lift(key)?.split(',').map((s) => s.trim()).filter(Boolean);
    const m: Row = {};
    switch (type) {
      case 'app':
        Object.assign(m, {browse_path: lift('browsePath'), icon: lift('icon'), url: lift('url'), admin: roles.includes('adminApp') ? true : undefined});
        break;
      case 'panel':
        m.condition = header.keys.condition?.[0];
        if (m.condition !== undefined) lifted.add('condition');
        m.target_type = header.inputs[0]?.type;
        if (header.inputs[0]?.options.semType) m.target_semtype = semtypeId(header.inputs[0].options.semType);
        break;
      case 'viewer':
        Object.assign(m, {icon: lift('icon'), trellisable: flag('trellisable'), grid_chart: flag('gridChart')});
        break;
      case 'filter':
        if (meta.semType) {
          m.semtype = semtypeId(lift('semType'));
          uses.push({role: 'filters', semtype: meta.semType, derived_by: 'annotation', confidence: 1});
        }
        Object.assign(m, {primary: flag('primaryFilter'), columnless: flag('columnlessFilter')});
        break;
      case 'cell-renderer':
        if (!meta.cellType) return null;
        m.cell_type = lift('cellType');
        m.column_tags = list('columnTags');
        uses.push({role: 'renders', semtype: meta.cellType, derived_by: 'annotation', confidence: 1});
        break;
      case 'sem-type-detector':
        if (meta.skipTest !== undefined) {
          lift('skipTest');
          m.skip_test = true;
        }
        if (meta.semType) uses.push({role: 'detects', semtype: lift('semType'), derived_by: 'annotation', confidence: 1});
        break;
      case 'file-handler':
        m.extensions = list('ext');
        m.direction = roles.includes('fileHandler') || roles.includes('file-handler') ? 'import' : 'export';
        break;
      case 'file-viewer':
        if (!meta.fileViewer) return null;
        m.extensions = list('fileViewer');
        m.check = lift('fileViewerCheck');
        break;
      case 'script-handler':
        if (!meta['scriptHandler.language'] || !meta['scriptHandler.extensions']) return null;
        Object.assign(m, {script_language: lift('scriptHandler.language'), extensions: list('scriptHandler.extensions'),
          comment_start: lift('scriptHandler.commentStart'), editor_mode: lift('scriptHandler.codeEditorMode')});
        break;
      case 'lifecycle-hook':
        m.phase = roles.includes('init') ? 'init' : 'autostart';
        m.immediate = flag('autostartImmediate');
        break;
    }
    return {members: m, consumed, uses};
  }

  private emitFunction(row: Row, uses: SemtypeUse[], header: Header, pkg: PackageFolder, file: string): void {
    const id = row.id as string;
    this.emitter.node(row);
    this.declares(pkg, id, file);
    for (const u of uses) {
      this.emitter.node({type: 'semantic-type', id: semtypeId(u.semtype), name: u.semtype, language: 'other', provenance: u.derived_by, source_layer: 'public'});
      this.emitter.edge({type: 'targets-semtype', from: id, to: semtypeId(u.semtype), role: u.role, derived_by: u.derived_by, confidence: u.confidence, evidence: [file]});
    }
    const feature = (header.keys.feature?.[0] ?? header.meta.feature)?.trim().replace(/^~/, '');
    if (!feature) return;
    this.emitter.edge({type: 'is-implemented-in', from: feature, to: id, derived_by: 'annotation', confidence: 1, evidence: [file]});
    this.emitter.claim({file, feature, rung: 1, source: 'marker', props: {}, line: header.line});
  }

  /** `meta.semType` failing, what a detector returns or assigns: `DG.SEMTYPE.X`, a string literal, or a file-level constant. */
  private detectedSemtypes(block: HeaderBlock, text: string, file: string): string[] {
    const out = new Set<string>();
    for (const pattern of DETECTED)
      for (const m of block.body.matchAll(pattern)) {
        if (m[1] && this.semtypes[m[1]]) out.add(this.semtypes[m[1]]);
        else if (m[3]) out.add(m[3]);
        else if (m[4]) {
          const value = constantValue(text, m[4], m[5]);
          if (value) out.add(value);
          else this.emitter.problem('unresolved_ids', `${file}: ${m[4]}${m[5] ? `.${m[5]}` : ''} names no semantic type this file declares`);
        }
      }
    return [...out];
  }

  private fileNode(pkg: PackageFolder, file: string): void {
    if (this.files.has(file)) return;
    this.files.add(file);
    this.emitter.node({type: 'source-file', id: fileId(file), name: path.posix.basename(file), path: file, loc: countLines(path.join(this.repoRoot, file)), language: languageOf(file),
      generated: file.endsWith('.g.ts') ? true : undefined, package: pkgId(pkg.folder), provenance: 'filesystem', source_layer: 'public'});
    this.declares(pkg, fileId(file), file);
  }

  private declares(pkg: PackageFolder, to: string, evidence: string): void {
    this.emitter.edge({type: 'declares', from: pkgId(pkg.folder), to, derived_by: 'registry', confidence: 1, evidence: [evidence]});
  }

  /** A name answers to itself and, when nothing matches it exactly, to its spaceless form. */
  private register(pkg: string, name: string, id: string): void {
    add(this.exact, exactKey(pkg, name), id);
    add(this.loose, looseKey(pkg, name), id);
  }

  private glob(pkg: PackageFolder, pattern: string): string[] {
    return globSync(`${pkg.dir}/${pattern}`, {cwd: this.repoRoot, ignore: SOURCE_IGNORE, nodir: true, posix: true, windowsPathsNoEscape: true}).sort();
  }

  private exists(file: string): boolean {
    return fs.existsSync(path.join(this.repoRoot, file));
  }

  private read(file: string): string {
    return fs.readFileSync(path.join(this.repoRoot, file), 'utf8');
  }
}

/** The package matches case-insensitively, the name too; `Pareto Front` and `Pareto front` share this key, the declaration
 * names `paretoFront` and `paretoFrontViewer` do not. */
function exactKey(pkg: string, name: string): string {
  return `${pkg.toLowerCase()}:${name.toLowerCase()}`;
}

/** `EDA:ParetoFront` as the generated package-api writes it, matched against the same name without its spaces. */
function looseKey(pkg: string, name: string): string {
  return exactKey(pkg, name).replace(/\s+/g, '');
}

function add(index: Map<string, Set<string>>, key: string, id: string): void {
  let ids = index.get(key);
  if (!ids) index.set(key, ids = new Set());
  ids.add(id);
}

/** `SEMTYPEGIS.GISCOUNTRY` or a plain `MOLECULE`: the string the constant holds, when its literal is in the same file. */
function constantValue(text: string, name: string, key?: string): string | undefined {
  if (!key) return new RegExp(`\\b${name}\\s*=\\s*(['"])([^'"\\n]+)\\1`).exec(text)?.[2];
  const literal = new RegExp(`\\b${name}\\s*=\\s*\\{`).exec(text);
  if (!literal) return undefined;
  const body = text.slice(literal.index + literal[0].length);
  const end = body.indexOf('}');
  return new RegExp(`\\b${key}\\s*:\\s*(['"])([^'"\\n]+)\\1`).exec(end < 0 ? body : body.slice(0, end))?.[2];
}

/** Conda dependency names, the nested `{pip: [...]}` lists included, without version pins. */
function dependencyNames(deps: unknown): string[] {
  const out: string[] = [];
  const visit = (item: unknown) => {
    if (typeof item === 'string') out.push(item.split(/[=<>]/)[0].trim());
    else if (Array.isArray(item)) item.forEach(visit);
    else if (item && typeof item === 'object') Object.values(item).forEach(visit);
  };
  visit(deps);
  return [...new Set(out.filter(Boolean))];
}

/** `DG.SEMTYPE.MOLECULE` -> `Molecule`, from js-api/src/const.ts; empty when the file is absent. */
function loadSemtypes(repoRoot: string): Record<string, string> {
  const file = path.join(repoRoot, 'public', 'js-api', 'src', 'const.ts');
  if (!fs.existsSync(file)) return {};
  const text = fs.readFileSync(file, 'utf8');
  const start = text.indexOf('export const SEMTYPE');
  if (start < 0) return {};
  const block = text.slice(start, text.indexOf('};', start));
  return Object.fromEntries([...block.matchAll(/^\s*([A-Za-z][A-Za-z0-9_]*):\s*'([^']*)'/gm)].map((m) => [m[1], m[2]]));
}
