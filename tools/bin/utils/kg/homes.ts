/// Home documents: the markdown files and YAML records whose frontmatter declares a node
/// (conventions.md §5), and the pages that only annotate one (`documents:`). Discovers them
/// under the monorepo, validates them against the type system, resolves references and
/// citations, and reports what `grok kg gen` would stub.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {splitFrontmatter, parseYamlDocument, keyLine, Frontmatter} from './frontmatter';
import {extractCitations, headings, proseLines, Citation} from './citations';
import {TypeSystem, NodeType, EdgeType, Member, Issue, checkValue, isSubtype, pascal, kebabOfLabel, concreteAuthored} from './types';
import {normalizeRow} from './normalize';
import {SCHEME_TYPES, PREFIXED_ID, SCHEMED_ID, PAGE_PATH, docCandidates} from './ids';

/** Where homes may live, relative to the monorepo root: any markdown document in the repos, and
 * the YAML records (concepts, people, teams, customers) inside the knowledge-graph folder. */
export const HOME_ROOTS = [
  'core/**/*.{md,mdx}',
  'public/**/*.{md,mdx}',
  'infra/**/*.{md,mdx}',
  'landing/**/*.{md,mdx}',
  'core/docs/knowledge-graph/**/*.yaml',
];

export const HOME_IGNORE = [
  'core/docs/knowledge-graph/nodes/**',
  'core/docs/knowledge-graph/edges/**',
  'core/docs/knowledge-graph/schema.yaml',
  'core/docs/knowledge-graph/questions/**',
  '**/node_modules/**',
  '**/dist/**',
  '**/.dart_tool/**',
  '**/build/**',
  '**/.claude/**',
  '**/.git/**',
  // test fixtures carry frontmatter of their own (grok-core mail/report fixtures, the kg fixtures)
  '**/fixtures/**',
  '**/__tests__/**',
];

/** Test Track scenarios (and their Playwright copies): only files migrated to `id: TS:...` are homes; the
 * legacy `feature:` key there still means the area (nodes/scenario.yaml rules). */
const LEGACY_TEST_TRACK = ['public/packages/UsageAnalysis/files/', 'public/playwright-public/'];
const INTERNAL_DIR = 'core/docs/knowledge-graph/internal/';
const DOCUSAURUS_KEYS = ['title', 'description', 'keywords', 'sidebar_position', 'sidebar_label', 'slug', 'mdx',
  'unlisted', 'toc_max_heading_level', 'position', 'format', 'pagination_prev', 'pagination_next', 'hide_sidebar',
  'hide_search', 'hide_title', 'tags'];
const BARE_ID = /^[a-z0-9-]+(\/[a-z0-9-]+)*$/;
const SEGMENT = /^[a-z0-9-]+$/;
export const REPO_PREFIX = /^(landing|infra):(.+)$/;
export const GLOB_MAGIC = /[*?[\]{}]/;

export interface Home {
  /** With the type prefix, without the sigil. */
  id: string;
  prefix?: string;
  /** Without the prefix. */
  localId: string;
  type: NodeType;
  /** Posix path relative to the monorepo root. */
  file: string;
  /** A YAML record rather than a markdown document: no body, `name:` required. */
  yaml: boolean;
  line: number;
  name: string;
  aliases: string[];
  data: Record<string, unknown>;
  body: string;
  fm: Frontmatter;
}

/** A page that is not a home but carries edge keys such as `documents:`. */
export interface AnnotatedPage {
  file: string;
  fm: Frontmatter;
}

export interface Stub {
  id: string;
  type: NodeType;
  neededBy: string;
}

/** A reference whose target is an extracted node the build has not produced yet. */
export interface UnresolvedRef {
  source: string;
  key: string;
  value: string;
  expectedTypes: string[];
}

export interface HomeSet {
  homes: Home[];
  pages: AnnotatedPage[];
  stubs: Stub[];
  errors: Issue[];
  warnings: Issue[];
  unresolvedExternal: UnresolvedRef[];
  /** Files walked looking for homes. */
  scanned: number;
  /** Pages that are not homes but carry an edge key such as `documents:`. */
  annotatedPages: number;
  citations: {doc: number, code: number, media: number};
  /** The homes by id and by alias, built once here and shared by every reader. */
  index: HomeIndex;
}

/** The `grok kg check` diagnostic report (not the graph's report node). */
export interface CheckReport {
  errors: Issue[];
  warnings: Issue[];
  /** Issue count per code. */
  codes: Record<string, number>;
  types: {nodes: number, edges: number, prefixes: number};
  scanned: number;
  homes: Record<string, number>;
  annotatedPages: number;
  citations: {doc: number, code: number, media: number};
  unresolvedExternal: UnresolvedRef[];
  stubs: string[];
  stale?: string[];
  written?: string[];
}

export function discoverHomeFiles(repoRoot: string): string[] {
  return globSync(HOME_ROOTS, {cwd: repoRoot, ignore: HOME_IGNORE, nodir: true, posix: true}).sort();
}

export function loadHomes(system: TypeSystem, repoRoot: string, files: string[] = discoverHomeFiles(repoRoot)): HomeSet {
  const set: HomeSet = {homes: [], pages: [], stubs: [], errors: [], warnings: [], unresolvedExternal: [], scanned: files.length,
    annotatedPages: 0, citations: {doc: 0, code: 0, media: 0}, index: {byId: new Map(), byAlias: new Map()}};
  for (const file of files) {
    let text: string;
    try {
      text = fs.readFileSync(path.join(repoRoot, file), 'utf8');
    } catch (e: any) {
      set.errors.push({file, code: 'unreadable', message: `cannot read: ${e.message}`});
      continue;
    }
    const isYaml = /\.ya?ml$/.test(file);
    const fm = isYaml ? parseYamlDocument(text) : splitFrontmatter(text);
    const legacy = LEGACY_TEST_TRACK.some((p) => file.startsWith(p));
    if (fm.error) {
      if (!legacy || /^id:\s*TS:/m.test(fm.yaml))
        set.errors.push({file, line: fm.errorLine, code: /not closed/.test(fm.error) ? 'unclosed-frontmatter' : 'yaml-error', message: fm.error});
      continue;
    }
    if (!fm.data) continue;
    if (legacy && !(typeof fm.data.id === 'string' && fm.data.id.startsWith('TS:'))) continue;
    if (fm.data.feature === undefined && fm.data.id === undefined) {
      if (isYaml) set.warnings.push({file, line: 1, code: 'stray-yaml', message: 'stray YAML file in the knowledge-graph folder: no id: or feature: key, so not a home'});
      else if (Object.keys(fm.data).some((k) => k === 'edges' || system.keys.has(k))) set.pages.push({file, fm});
      continue;
    }
    const home = readHome(system, file, fm, isYaml, set.errors);
    if (home) set.homes.push(home);
  }
  set.index = indexHomes(set.homes, set.errors);
  const checker = new HomeChecker(system, repoRoot, set.index, set);
  for (const home of set.homes)
    checker.check(home, home.fm);
  for (const page of set.pages)
    checker.checkPage(page.file, page.fm);
  checker.checkCycles();
  return set;
}

function readHome(system: TypeSystem, file: string, fm: Frontmatter, isYaml: boolean, errors: Issue[]): Home | null {
  const data = fm.data!;
  const error = (code: string, message: string, key = 'id') => errors.push({file, line: keyLine(fm, key) ?? keyLine(fm, 'feature') ?? 1, code, message});
  if (data.feature !== undefined && data.id !== undefined) {
    error('bad-id', 'both feature: and id: are set; feature: is sugar for id: without a prefix');
    return null;
  }
  const key = data.feature !== undefined ? 'feature' : 'id';
  const raw = data[key];
  if (typeof raw !== 'string' || !raw.trim()) {
    error('bad-id', `${key}: must be an id string, got ${JSON.stringify(raw)}`, key);
    return null;
  }
  const id = raw.trim().replace(/^~/, '');
  const prefixed = PREFIXED_ID.exec(id);
  if (key === 'feature' && prefixed) {
    error('bad-id', `feature: takes a feature id without a prefix; use id: ${id}`, key);
    return null;
  }
  const prefix = prefixed ? prefixed[1] : undefined;
  const localId = prefixed ? prefixed[2] : id;
  const baseName = prefix ? system.prefixes.get(prefix) : 'feature';
  if (!baseName || !system.nodes.has(baseName)) {
    error('bad-id', `unknown id prefix '${prefix}'`, key);
    return null;
  }
  let type = system.nodes.get(baseName)!;
  if (data.type !== undefined) {
    const named = system.nodes.get(String(data.type));
    if (!named) {
      error('bad-type', `unknown type '${data.type}'`, 'type');
      return null;
    }
    if (!isSubtype(system, named.name, baseName)) {
      error('bad-type', `type '${named.name}' is not a subtype of ${baseName}, the type of prefix ${prefix ?? '(none)'}`, 'type');
      return null;
    }
    type = named;
  }
  if (type.abstract) {
    error('bad-type', `type '${type.name}' is abstract; set type: to a concrete subtype`, data.type !== undefined ? 'type' : key);
    return null;
  }
  if (!type.authored) {
    error('bad-type', `type '${type.name}' is extracted, not authored; it cannot have a home document`, key);
    return null;
  }
  const segments = localId.split('/');
  if (!segments.every((s) => SEGMENT.test(s))) {
    error('bad-id', `id '${id}': segments must be lowercase kebab-case separated by /`, key);
    return null;
  }
  if (segments.length > 1 && !type.hierarchical) {
    error('bad-id', `id '${id}': type ${type.name} is not hierarchical, the id may not contain /`, key);
    return null;
  }
  const aliases = Array.isArray(data.aliases) ? data.aliases.map((a) => {
    const alias = String(a).replace(/^~/, '');
    return prefix && !PREFIXED_ID.test(alias) ? `${prefix}:${alias}` : alias;
  }) : [];
  const name = typeof data.name === 'string' ? data.name :
    isYaml ? undefined : typeof data.title === 'string' ? data.title : firstHeading(fm.body);
  return {id, prefix, localId, type, file, yaml: isYaml, line: keyLine(fm, key) ?? 1, name: name ?? '', aliases, data, body: fm.body, fm};
}

export function firstHeading(body: string): string | undefined {
  for (const {text} of proseLines(body)) {
    const m = /^#\s+(.+?)\s*#*\s*$/.exec(text);
    if (m) return m[1];
  }
  return undefined;
}

export interface HomeIndex {
  byId: Map<string, Home>;
  byAlias: Map<string, Home>;
}

/** The first of [candidates] that names a home, by id or by alias. */
export function lookupHome(index: HomeIndex, candidates: string[]): Home | undefined {
  for (const c of candidates) {
    const home = index.byId.get(c) ?? index.byAlias.get(c);
    if (home) return home;
  }
  return undefined;
}

function indexHomes(homes: Home[], errors: Issue[]): HomeIndex {
  const byId = new Map<string, Home>();
  const byAlias = new Map<string, Home>();
  for (const home of homes) {
    const other = byId.get(home.id);
    if (!other) {
      byId.set(home.id, home);
      continue;
    }
    // discovery order says nothing about which file is the newer claim, so both are named
    errors.push({file: home.file, line: home.line, code: 'duplicate-id', message: `duplicate id ~${home.id}, also the home of ${other.file}`});
    errors.push({file: other.file, line: other.line, code: 'duplicate-id', message: `duplicate id ~${home.id}, also the home of ${home.file}`});
  }
  for (const home of homes)
    for (const alias of home.aliases) {
      const owner = byId.get(alias);
      if (owner) {
        errors.push({file: home.file, line: home.line, code: 'alias-conflict', message: `alias ~${alias} is the id of ${owner.file}`});
        continue;
      }
      const other = byAlias.get(alias);
      if (other && other !== home) errors.push({file: home.file, line: home.line, code: 'alias-conflict', message: `alias ~${alias} is also an alias of ${other.file}`});
      else byAlias.set(alias, home);
    }
  return {byId, byAlias};
}

/** What an edge key hangs off: a home, or a plain page (a doc-page) that only annotates. */
interface Subject {
  file: string;
  typeName: string;
  home?: Home;
}

type ErrorFn = (code: string, message: string, key?: string, target?: string) => void;

class HomeChecker {
  private pathCache = new Map<string, string | null>();
  private anchorCache = new Map<string, Set<string> | null>();
  private stubIds = new Set<string>();
  /** Instances of `acyclic` edges among homes, edge -> from id -> to ids. */
  private instances = new Map<string, Map<string, Set<string>>>();

  constructor(private system: TypeSystem, private repoRoot: string, private index: HomeIndex, private set: HomeSet) {}

  check(home: Home, fm: Frontmatter): void {
    const error: ErrorFn = (code, message, key, target) => this.set.errors.push({file: home.file, line: key ? keyLine(fm, key) ?? home.line : home.line, code, message, target});
    const warn: ErrorFn = (code, message, key) => this.set.warnings.push({file: home.file, line: key ? keyLine(fm, key) ?? home.line : home.line, code, message});
    const {type, data} = home;
    const subject: Subject = {file: home.file, typeName: type.name, home};
    const isHelp = home.file.startsWith('public/help/');
    if (!home.name && (data.name === undefined || data.name === null))
      error('no-name', home.yaml ? 'no name: a YAML home must set name:' : 'no name: set name:, title:, or start the body with a # heading');
    const row: Record<string, unknown> = {type: type.name};
    for (const [key, value] of Object.entries(data)) {
      if (key === 'feature' || key === 'id' || key === 'type' || value === null) continue;
      if (key === 'title' && typeof value === 'string') continue;
      if (key === 'part_of') {
        error('part-of-authored', 'part-of is derived from the id path and never authored; remove part_of:', key);
        continue;
      }
      if (key === 'edges') {
        this.checkEdgesKey(subject, value, error);
        continue;
      }
      const edge = this.system.keys.get(key);
      if (edge) {
        this.checkEdgeKey(subject, key, edge, value, error, warn);
        continue;
      }
      const member = type.members[key];
      if (member) {
        if (isHomeFileMember(member) && !home.yaml) error('authored-path', `${key}: derived from the home file (${home.file}) and never authored; remove it`, key);
        else row[key] = value;
        continue;
      }
      if (isHelp && DOCUSAURUS_KEYS.includes(key)) continue;
      error('unknown-key', `unknown key '${key}' for type ${type.name}`, key);
    }
    const normalized = normalizeRow(this.system, row, {
      path: (v) => this.pathProblem(v),
      ref: (v, m) => this.resolveRef(v, m.refs!, {source: home.file, key: m.name}).problem ?? null,
    });
    for (const problem of normalized.problems) error(problem.code, problem.message, problem.key, problem.target);
    for (const member of Object.values(type.members)) {
      if (member.nullable || member.name === 'id' || member.name === 'name') continue;
      if (isHomeFileMember(member) && !home.yaml) {
        data[member.name] = home.file;
        continue;
      }
      if (normalized.row[member.name] === undefined) error('missing-key', `missing required key '${member.name}' for type ${type.name}`);
    }
    if (type.members.manual_only && data.manual_only === undefined && data.target_layer === 'manual-only')
      data.manual_only = true;
    if (typeof data.visibility === 'string' && home.file.startsWith(INTERNAL_DIR) && data.visibility !== 'internal')
      error('visibility-ceiling', `visibility '${data.visibility}' exceeds the ceiling 'internal' of ${INTERNAL_DIR}; a home there can only be internal`, 'visibility');
    this.checkHierarchy(home, error);
    if (!home.yaml) this.checkCitations(home, fm);
  }

  /** A page that is not a home but carries edge keys: `documents:` on a help page. */
  checkPage(file: string, fm: Frontmatter): void {
    const error: ErrorFn = (code, message, key, target) => this.set.errors.push({file, line: key ? keyLine(fm, key) ?? 1 : 1, code, message, target});
    const warn: ErrorFn = (code, message, key) => this.set.warnings.push({file, line: key ? keyLine(fm, key) ?? 1 : 1, code, message});
    if (!this.system.nodes.has('doc-page')) {
      error('bad-edge', 'a page can only annotate when the type system has a doc-page type');
      return;
    }
    const subject: Subject = {file, typeName: 'doc-page'};
    this.set.annotatedPages++;
    for (const [key, value] of Object.entries(fm.data!)) {
      if (value === null) continue;
      if (key === 'edges') this.checkEdgesKey(subject, value, error);
      else {
        const edge = this.system.keys.get(key);
        if (edge) this.checkEdgeKey(subject, key, edge, value, error, warn);
      }
    }
  }

  private checkHierarchy(home: Home, error: ErrorFn): void {
    const {type, localId} = home;
    if (!type.hierarchical) return;
    const segments = localId.split('/');
    const withPrefix = (n: number) => (home.prefix ? `${home.prefix}:` : '') + segments.slice(0, n).join('/');
    const missing = (n: number) => !this.index.byId.has(withPrefix(n));
    let areaMissing = false;
    if (type.root === 'feature') {
      const roots = this.system.featureRoots;
      if (roots.length && !roots.includes(segments[0]))
        error('bad-root', `~${home.id}: feature root '${segments[0]}' is not one of ${roots.join(', ')} (schema.yaml feature_roots)`);
      if (segments.length === 2 && (home.data.owner === undefined || home.data.owner === null))
        error('missing-owner', `~${home.id} is a second-level feature node and must set owner:`);
      if (segments.length >= 3 && missing(2)) {
        areaMissing = true;
        error('missing-area', `~${home.id} needs its area home ~${withPrefix(2)}: a level-2 feature must exist and carry the owner; none found`);
      }
    }
    if (type.root === 'concept' && segments.length > 2)
      error('too-deep', `~${home.id}: a concept id has at most two segments`);
    for (let i = 1; i < segments.length; i++) {
      const ancestor = withPrefix(i);
      if (!missing(i) || this.stubIds.has(ancestor) || (areaMissing && i === 2)) continue;
      this.stubIds.add(ancestor);
      this.set.stubs.push({id: ancestor, type, neededBy: home.id});
      this.set.warnings.push({file: home.file, line: home.line, code: 'stub', message: `no home for ~${ancestor} (parent of ~${home.id}); grok kg gen lists it as a stub`});
    }
  }

  private checkEdgeKey(subject: Subject, key: string, edge: EdgeType, value: unknown, error: ErrorFn, warn: ErrorFn): void {
    if (edge.abstract) {
      error('bad-edge', `${key}: spells the abstract edge ${edge.name}; abstract edges cannot be authored`, key);
      return;
    }
    if (!Array.isArray(value)) {
      error('bad-edge', `${key}: expected a list of ${edge.name} targets, got ${JSON.stringify(value)}`, key);
      return;
    }
    const homeSide = edge.keySide;
    const otherSide = homeSide === 'from' ? 'to' : 'from';
    if (!edge[homeSide].some((t) => isSubtype(this.system, subject.typeName, t))) {
      error('bad-edge', `${key}: ${edge.name}.${homeSide} must be ${edge[homeSide].map(pascal).join(' | ')}; this ${subject.home ? 'home' : 'page'} is a ${subject.typeName}`, key);
      return;
    }
    if (homeSide === 'from' && edge.cardinality === 'one' && value.length > 1)
      error('cardinality', `${key}: ${edge.name} has cardinality one, ${value.length} targets given`, key);
    const targetKey = key === 'code' ? 'path' : 'to';
    const seen = new Set<string>();
    value.forEach((item, i) => {
      const where = `${key}[${i}]`;
      let target: unknown;
      let props: Record<string, unknown> = {};
      if (typeof item === 'string') target = item;
      else if (item && typeof item === 'object' && !Array.isArray(item)) {
        ({[targetKey]: target, ...props} = item as Record<string, unknown>);
        if (target === undefined) {
          error('bad-edge', `${where}: a map item needs ${targetKey}:`, key);
          return;
        }
      } else {
        error('bad-edge', `${where}: expected an id or a map with ${targetKey}:, got ${JSON.stringify(item)}`, key);
        return;
      }
      if (typeof target !== 'string' || !target.trim()) {
        error('bad-edge', `${where}: ${targetKey} must be a string, got ${JSON.stringify(target)}`, key);
        return;
      }
      if (seen.has(target.trim())) warn('duplicate-item', `${where}: '${target.trim()}' is listed twice under ${key}:`, key);
      seen.add(target.trim());
      if (key === 'code') {
        const problem = this.pathProblem(target);
        if (problem) error('missing-path', `${where}: ${problem}`, key, target);
      } else {
        const {home, problem} = this.resolveRef(target, edge[otherSide], {source: subject.file, key: where});
        if (problem) error('unresolved-ref', `${where}: ${problem}`, key);
        if (home && subject.home) {
          if (edge.sameType && home.type.name !== subject.typeName)
            error('same-type', `${where}: ${edge.name} is same_type; ~${home.id} is a ${home.type.name}, this home is a ${subject.typeName}`, key);
          if (edge.acyclic) this.addInstance(edge, homeSide === 'from' ? subject.home.id : home.id, homeSide === 'from' ? home.id : subject.home.id);
        }
      }
      this.checkEdgeProperties(edge, props, where, error, key);
    });
  }

  private checkEdgesKey(subject: Subject, value: unknown, error: ErrorFn): void {
    if (!Array.isArray(value)) {
      error('bad-edge', `edges: expected a list of {type, to, ...} maps, got ${JSON.stringify(value)}`, 'edges');
      return;
    }
    value.forEach((item, i) => {
      const where = `edges[${i}]`;
      if (!item || typeof item !== 'object' || Array.isArray(item)) {
        error('bad-edge', `${where}: expected a map with type: and to:`, 'edges');
        return;
      }
      const {type, to, ...props} = item as Record<string, unknown>;
      const edge = typeof type === 'string' ? this.system.edges.get(type) : undefined;
      if (!edge) {
        const kebab = typeof type === 'string' ? kebabOfLabel(type) : null;
        if (kebab && this.system.edges.has(kebab)) error('bad-edge', `${where}: edge types are lower-dash-case: write '${kebab}', not '${type}'`, 'edges');
        else error('bad-edge', `${where}: unknown edge type ${JSON.stringify(type)}`, 'edges');
        return;
      }
      if (edge.abstract || !edge.derivedBy.includes('annotation')) {
        error('bad-edge', `${where}: ${edge.name} is ${edge.abstract ? 'abstract' : 'never authored (derived_by lacks annotation)'}`, 'edges');
        return;
      }
      if (!edge.from.some((t) => isSubtype(this.system, subject.typeName, t))) {
        error('bad-edge', `${where}: ${edge.name}.from must be ${edge.from.map(pascal).join(' | ')}; this ${subject.home ? 'home' : 'page'} is a ${subject.typeName}`, 'edges');
        return;
      }
      if (typeof to !== 'string' || !to.trim()) {
        error('bad-edge', `${where}: to must be an id, got ${JSON.stringify(to)}`, 'edges');
        return;
      }
      const {home, problem} = this.resolveRef(to, edge.to, {source: subject.file, key: where});
      if (problem) error('unresolved-ref', `${where}: ${problem}`, 'edges');
      if (home && subject.home) {
        if (edge.sameType && home.type.name !== subject.typeName)
          error('same-type', `${where}: ${edge.name} is same_type; ~${home.id} is a ${home.type.name}, this home is a ${subject.typeName}`, 'edges');
        if (edge.acyclic) this.addInstance(edge, subject.home.id, home.id);
      }
      this.checkEdgeProperties(edge, props, where, error, 'edges');
    });
  }

  private checkEdgeProperties(edge: EdgeType, props: Record<string, unknown>, where: string, error: ErrorFn, key: string): void {
    for (const [name, value] of Object.entries(props)) {
      const member = edge.properties[name];
      if (!member) {
        error('unknown-key', `${where}: ${edge.name} has no property '${name}' (${Object.keys(edge.properties).join(', ') || 'none'})`, key);
        continue;
      }
      if (value === null) continue;
      const problem = this.valueProblem(member, value, {source: '', key: `${where}.${name}`});
      if (problem) error('bad-value', `${where}.${name}: ${problem}`, key);
    }
    for (const member of Object.values(edge.properties))
      if (!member.nullable && (props[member.name] === undefined || props[member.name] === null))
        error('missing-key', `${where}: ${edge.name} requires '${member.name}'`, key);
  }

  private addInstance(edge: EdgeType, from: string, to: string): void {
    let byFrom = this.instances.get(edge.name);
    if (!byFrom) this.instances.set(edge.name, byFrom = new Map());
    let tos = byFrom.get(from);
    if (!tos) byFrom.set(from, tos = new Set());
    tos.add(to);
  }

  /** After every home is read: an `acyclic` edge whose instances form a cycle is an error on the home closing it. */
  checkCycles(): void {
    for (const [edgeName, byFrom] of this.instances) {
      const state = new Map<string, 'open' | 'done'>();
      const visit = (id: string, trail: string[]): void => {
        state.set(id, 'open');
        for (const to of byFrom.get(id) ?? []) {
          const s = state.get(to);
          if (s === 'done') continue;
          if (s === 'open') {
            const cycle = [...trail.slice(trail.indexOf(to)), id, to];
            const home = this.index.byId.get(id)!;
            this.set.errors.push({file: home.file, line: home.line, code: 'acyclic', message: `${edgeName} is acyclic but forms a cycle: ${cycle.map((c) => `~${c}`).join(' -> ')}`});
            continue;
          }
          visit(to, [...trail, id]);
        }
        state.set(id, 'done');
      };
      for (const id of byFrom.keys())
        if (!state.has(id)) visit(id, []);
    }
  }

  private valueProblem(member: Member, value: unknown, ctx: {source: string, key: string}): string | null {
    return checkValue(member, value, {
      provenance: this.system.provenance,
      path: (v) => this.pathProblem(v),
      ref: (v, m) => this.resolveRef(v, m.refs!, ctx).problem ?? null,
    });
  }

  /** A path or glob relative to the monorepo root, `landing:`/`infra:` prefixed for those repos, `#anchor` allowed. */
  pathProblem(value: string): string | null {
    const cached = this.pathCache.get(value);
    if (cached !== undefined) return cached;
    if (value.includes('\\')) return `path '${value}' contains a backslash; use forward slashes`;
    let p = value.split('#')[0].trim();
    const repo = REPO_PREFIX.exec(p);
    if (repo) p = `${repo[1]}/${repo[2]}`;
    p = p.replace(/\/+$/, '');
    let problem: string | null;
    if (GLOB_MAGIC.test(p))
      problem = globSync(p, {cwd: this.repoRoot, ignore: HOME_IGNORE, posix: true, windowsPathsNoEscape: true}).length ? null : `no file matches '${value}'`;
    else
      problem = fs.existsSync(path.join(this.repoRoot, p)) ? null : `path '${value}' does not exist`;
    this.pathCache.set(value, problem);
    return problem;
  }

  /**
   * Resolves a reference written in a home: parse the id, expand the expected types to the concrete
   * authored types that may carry it, check the id's kind against them, then look it up. Extracted ids
   * (schemes, tracker keys, extracted prefixes) are type-checked now and listed as unresolved-external.
   */
  resolveRef(value: string, expected: string[], ctx: {source: string, key: string}): {home?: Home, problem?: string} {
    const raw = value.trim().replace(/^~/, '');
    const wanted = expected.map(pascal).join(' | ');
    const compatible = (types: string[]) => types.some((t) => expected.some((e) => isSubtype(this.system, t, e)));
    const external = (types: string[]) => {
      if (!compatible(types)) return {problem: `'${value}' is a ${types.map(pascal).join(' | ')}; expected ${wanted}`};
      this.set.unresolvedExternal.push({source: ctx.source, key: ctx.key, value, expectedTypes: expected.map(pascal)});
      return {};
    };
    if (!raw) return {problem: 'empty id'};
    if (/^GROK-\d+$/.test(raw)) return external(['ticket'].filter((t) => this.system.nodes.has(t)));
    if (/^#\d+$/.test(raw)) return external(['ticket', 'pull-request'].filter((t) => this.system.nodes.has(t)));
    const prefixed = PREFIXED_ID.exec(raw);
    if (prefixed) {
      const typeName = this.system.prefixes.get(prefixed[1]);
      if (!typeName) return {problem: `'${value}': unknown id prefix '${prefixed[1]}'`};
      if (!compatible([typeName])) return {problem: `'${value}' is a ${typeName} (prefix ${prefixed[1]}); expected ${wanted}`};
      if (!this.system.nodes.get(typeName)!.authored) return external([typeName]);
      const {id, anchor} = splitRefAnchor(prefixed[2]);
      return this.lookup(value, [`${prefixed[1]}:${id}`], expected, anchor);
    }
    const schemed = SCHEMED_ID.exec(raw);
    if (schemed) {
      const types = (SCHEME_TYPES[schemed[1]] ?? []).filter((t) => this.system.nodes.has(t));
      if (!types.length) return {problem: `'${value}': unknown id scheme '${schemed[1]}:'`};
      if (schemed[1] === 'decl') {
        const problem = this.pathProblem(schemed[2].split('#')[0]);
        if (problem) return {problem: `'${value}': the declaration's ${problem}`};
      }
      return external(types);
    }
    if (this.system.nodes.has('doc-page') && PAGE_PATH.test(raw)) {
      const problem = this.pathProblem(raw);
      return problem ? {problem: `'${value}': the page's ${problem}`} : external(['doc-page']);
    }
    const {id, anchor} = splitRefAnchor(raw);
    if (!BARE_ID.test(id)) return {problem: `'${value}' is not a valid id: lowercase kebab segments, an optional Type: prefix, an optional #anchor`};
    const concrete = concreteAuthored(this.system, expected);
    if (!concrete.length) return {problem: `'${value}' cannot be a ${wanted}: an extracted id carries its scheme or tracker key (pkg:, func:, decl:, GROK-n, ...)`};
    const candidates = [...new Set(concrete.map((t) => t.prefix ? `${t.prefix}:${id}` : id))];
    return this.lookup(value, candidates, expected, anchor);
  }

  private lookup(value: string, candidates: string[], expected: string[], anchor?: string): {home?: Home, problem?: string} {
    const hits = [...new Set(candidates.map((c) => lookupHome(this.index, [c])).filter((h): h is Home => !!h))];
    if (hits.length > 1) return {problem: `'${value}' is ambiguous: ${hits.map((h) => `~${h.id}`).join(', ')}; write the prefix`};
    if (!hits.length) return {problem: `'${value}' does not resolve to any home document (as ${candidates.map((c) => `~${c}`).join(' or ')})`};
    const home = hits[0];
    if (!expected.some((t) => isSubtype(this.system, home.type.name, t)))
      return {problem: `'${value}' resolves to ~${home.id}, a ${home.type.name}; expected ${expected.map(pascal).join(' | ')}`};
    if (anchor !== undefined) {
      if (home.yaml) return {home, problem: `'${value}': ~${home.id} is a YAML record and has no headings, so '#${anchor}' cannot resolve`};
      if (!this.anchorsOf(home.file)?.has(anchor)) return {home, problem: `'${value}': no heading '#${anchor}' in ${home.file}`};
    }
    return {home};
  }

  /** Heading anchors of a markdown file under the repo root; null when it cannot be read. */
  private anchorsOf(file: string): Set<string> | null {
    const cached = this.anchorCache.get(file);
    if (cached !== undefined) return cached;
    let anchors: Set<string> | null;
    try {
      anchors = new Set(headings(splitFrontmatter(fs.readFileSync(path.join(this.repoRoot, file), 'utf8')).body).map((h) => h.slug));
    } catch {
      anchors = null;
    }
    this.anchorCache.set(file, anchors);
    return anchors;
  }

  private checkCitations(home: Home, fm: Frontmatter): void {
    for (const c of extractCitations(home.file, home.body, fm.bodyLine)) {
      const error = (code: string, message: string, target?: string) => this.set.errors.push({file: home.file, line: c.line, code, message, target});
      if (c.resolved === null) {
        this.set.citations.code++;
        error('citation-escape', `link '${c.raw}' escapes the repository`, c.raw);
        continue;
      }
      // a Docusaurus link may drop the extension: [Tile viewer](tile-viewer) means tile-viewer.md beside the page
      const extensionless = c.kind !== 'backtick' && !path.posix.extname(c.resolved) && this.pathProblem(c.resolved) !== null;
      const resolved = extensionless ? docCandidates(c.resolved).find((a) => !this.pathProblem(a)) ?? c.resolved : c.resolved;
      const target = resolved === c.resolved ? c.target : 'doc';
      this.set.citations[target]++;
      const problem = this.pathProblem(resolved);
      if (problem) {
        error(target === 'doc' ? 'missing-doc-link' : 'missing-cited-path', target === 'doc' ?
          `linked document '${c.raw}' does not exist${c.raw === resolved ? '' : ` (resolved to ${resolved})`}` : `cited ${problem}`,
          target === 'doc' ? c.raw : resolved);
        continue;
      }
      if (c.anchor && target === 'doc' && !this.anchorsOf(resolved)?.has(c.anchor))
        error('bad-anchor', `link '${c.raw}': no heading '#${c.anchor}' in ${resolved}`, c.raw);
    }
  }
}

/** A `path: Path` member on an authored type names the home file itself (scenario.path, doc-page.path). */
export function isHomeFileMember(member: Member): boolean {
  return member.name === 'path' && member.kind === 'scalar' && member.scalar === 'Path' && !member.list;
}

/** A reference's `#anchor` part; `#` inside an extracted id (`decl:file#Name`) is not an anchor. */
export function splitRefAnchor(id: string): {id: string, anchor?: string} {
  const hash = id.indexOf('#');
  return hash < 0 ? {id} : {id: id.slice(0, hash), anchor: id.slice(hash + 1)};
}

export function makeReport(system: TypeSystem, homes: HomeSet | null): CheckReport {
  const byType: Record<string, number> = {};
  for (const home of homes?.homes ?? [])
    byType[home.type.name] = (byType[home.type.name] ?? 0) + 1;
  const errors = [...system.errors, ...(homes?.errors ?? [])];
  const warnings = [...system.warnings, ...(homes?.warnings ?? [])];
  const codes: Record<string, number> = {};
  for (const issue of [...errors, ...warnings]) codes[issue.code] = (codes[issue.code] ?? 0) + 1;
  return {
    errors, warnings, codes,
    types: {nodes: system.nodes.size, edges: system.edges.size, prefixes: system.prefixes.size},
    scanned: homes?.scanned ?? 0,
    homes: byType,
    annotatedPages: homes?.annotatedPages ?? 0,
    citations: homes?.citations ?? {doc: 0, code: 0, media: 0},
    unresolvedExternal: homes?.unresolvedExternal ?? [],
    stubs: (homes?.stubs ?? []).map((s) => s.id),
  };
}

export type {Citation};
