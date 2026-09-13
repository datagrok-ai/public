/// Home documents: the markdown files whose frontmatter declares a node (CONVENTIONS §5).
/// Discovers them under the monorepo, validates the frontmatter against the type system,
/// resolves references and cited paths, and reports what `grok kg gen` would stub.
import * as fs from 'fs';
import * as path from 'path';
import {globSync} from 'glob';
import {splitFrontmatter, keyLine, Frontmatter} from './frontmatter';
import {TypeSystem, NodeType, EdgeType, Member, Issue, checkValue, isSubtype, pascal} from './types';

/** Where home documents may live, relative to the monorepo root: any markdown file in the repos. */
export const HOME_ROOTS = [
  'core/**/*.{md,mdx}',
  'public/**/*.{md,mdx}',
  'infra/**/*.{md,mdx}',
  'landing/**/*.{md,mdx}',
];

export const HOME_IGNORE = [
  '**/node_modules/**',
  '**/dist/**',
  '**/.dart_tool/**',
  '**/build/**',
  '**/.claude/**',
  '**/.git/**',
  // test fixtures carry frontmatter of their own (grok-core mail/report fixtures, the kg fixtures)
  '**/fixtures/**',
  '**/__tests__/**',
  // Test Track scenarios (and their Playwright copies) carry a `feature:` key with a different
  // meaning (the area they belong to) until they migrate to `covers:` (nodes/scenario.yaml rules).
  'public/packages/UsageAnalysis/files/**',
  'public/playwright-public/**',
];

const DOCUSAURUS_KEYS = ['title', 'description', 'keywords', 'sidebar_position', 'sidebar_label', 'slug', 'mdx',
  'unlisted', 'toc_max_heading_level', 'position', 'format', 'pagination_prev', 'pagination_next', 'hide_sidebar',
  'hide_search', 'hide_title', 'tags'];
const PREFIXED_ID = /^([A-Z][A-Za-z]{0,5}):(.+)$/;
const SEGMENT = /^[a-z0-9-]+$/;
const EXTERNAL_ID = /^([a-z][a-z0-9-]*:|GROK-\d+$|#\d+$)/;
const CITED_PATH = /^(core|public|infra|landing)[\\/][A-Za-z0-9_.\\/*-]+$/;
const REPO_PREFIX = /^(landing|infra):(.+)$/;
const GLOB_MAGIC = /[*?[\]{}]/;

export interface Home {
  /** With the type prefix, without the sigil. */
  id: string;
  prefix?: string;
  /** Without the prefix. */
  localId: string;
  type: NodeType;
  /** Posix path relative to the monorepo root. */
  file: string;
  line: number;
  name: string;
  aliases: string[];
  data: Record<string, unknown>;
  body: string;
}

export interface Stub {
  id: string;
  type: NodeType;
  neededBy: string;
}

export interface HomeSet {
  homes: Home[];
  stubs: Stub[];
  errors: Issue[];
  warnings: Issue[];
  unresolvedExternal: number;
  /** Markdown files walked looking for homes. */
  scanned: number;
}

export interface Report {
  errors: Issue[];
  warnings: Issue[];
  types: {nodes: number, edges: number, prefixes: number};
  scanned: number;
  homes: Record<string, number>;
  unresolvedExternal: number;
  stubs: string[];
  stale?: string[];
  written?: string[];
}

export function discoverHomeFiles(repoRoot: string): string[] {
  return globSync(HOME_ROOTS, {cwd: repoRoot, ignore: HOME_IGNORE, nodir: true, posix: true}).sort();
}

export function loadHomes(system: TypeSystem, repoRoot: string, files: string[] = discoverHomeFiles(repoRoot)): HomeSet {
  const set: HomeSet = {homes: [], stubs: [], errors: [], warnings: [], unresolvedExternal: 0, scanned: files.length};
  const fms = new Map<Home, Frontmatter>();
  for (const file of files) {
    const text = fs.readFileSync(path.join(repoRoot, file), 'utf8');
    const fm = splitFrontmatter(text);
    if (fm.error) {
      if (/^(feature|id):/m.test(fm.yaml)) set.errors.push({file, line: fm.errorLine, message: fm.error});
      continue;
    }
    if (!fm.data || (fm.data.feature === undefined && fm.data.id === undefined)) continue;
    const home = readHome(system, file, fm, set.errors);
    if (home) {
      set.homes.push(home);
      fms.set(home, fm);
    }
  }
  const index = indexHomes(set.homes, set.errors);
  const checker = new HomeChecker(system, repoRoot, index, set);
  for (const home of set.homes)
    checker.check(home, fms.get(home)!);
  return set;
}

function readHome(system: TypeSystem, file: string, fm: Frontmatter, errors: Issue[]): Home | null {
  const data = fm.data!;
  const error = (message: string, key = 'id') => errors.push({file, line: keyLine(fm, key) ?? keyLine(fm, 'feature') ?? 1, message});
  if (data.feature !== undefined && data.id !== undefined) {
    error('both feature: and id: are set; feature: is sugar for id: without a prefix');
    return null;
  }
  const key = data.feature !== undefined ? 'feature' : 'id';
  const raw = data[key];
  if (typeof raw !== 'string' || !raw.trim()) {
    error(`${key}: must be an id string, got ${JSON.stringify(raw)}`, key);
    return null;
  }
  const id = raw.trim().replace(/^~/, '');
  const prefixed = PREFIXED_ID.exec(id);
  if (key === 'feature' && prefixed) {
    error(`feature: takes a feature id without a prefix; use id: ${id}`, key);
    return null;
  }
  const prefix = prefixed ? prefixed[1] : undefined;
  const localId = prefixed ? prefixed[2] : id;
  const baseName = prefix ? system.prefixes.get(prefix) : 'feature';
  if (!baseName || !system.nodes.has(baseName)) {
    error(`unknown id prefix '${prefix}'`, key);
    return null;
  }
  let type = system.nodes.get(baseName)!;
  if (data.type !== undefined) {
    const named = system.nodes.get(String(data.type));
    if (!named) {
      error(`unknown type '${data.type}'`, 'type');
      return null;
    }
    if (!isSubtype(system, named.name, baseName)) {
      error(`type '${named.name}' is not a subtype of ${baseName}, the type of prefix ${prefix ?? '(none)'}`, 'type');
      return null;
    }
    type = named;
  }
  if (type.abstract) {
    error(`type '${type.name}' is abstract; set type: to a concrete subtype`, data.type !== undefined ? 'type' : key);
    return null;
  }
  if (!type.authored) {
    error(`type '${type.name}' is extracted, not authored; it cannot have a home document`, key);
    return null;
  }
  const segments = localId.split('/');
  if (!segments.every((s) => SEGMENT.test(s))) {
    error(`id '${id}': segments must be lowercase kebab-case separated by /`, key);
    return null;
  }
  if (segments.length > 1 && !type.hierarchical) {
    error(`id '${id}': type ${type.name} is not hierarchical, the id may not contain /`, key);
    return null;
  }
  const aliases = Array.isArray(data.aliases) ? data.aliases.map((a) => {
    const alias = String(a).replace(/^~/, '');
    return prefix && !PREFIXED_ID.test(alias) ? `${prefix}:${alias}` : alias;
  }) : [];
  const name = typeof data.name === 'string' ? data.name : typeof data.title === 'string' ? data.title : firstHeading(fm.body);
  return {id, prefix, localId, type, file, line: keyLine(fm, key) ?? 1, name: name ?? '', aliases, data, body: fm.body};
}

function firstHeading(body: string): string | undefined {
  for (const line of proseLines(body)) {
    const m = /^#\s+(.+?)\s*#*\s*$/.exec(line.text);
    if (m) return m[1];
  }
  return undefined;
}

/** Body lines outside fenced code blocks, with 0-based offsets. */
function proseLines(body: string): {text: string, offset: number}[] {
  const out: {text: string, offset: number}[] = [];
  let fence: string | null = null;
  body.split('\n').forEach((text, offset) => {
    const m = /^\s*(```|~~~)/.exec(text);
    if (m) {
      if (fence === null) fence = m[1];
      else if (fence === m[1]) fence = null;
      return;
    }
    if (fence === null) out.push({text, offset});
  });
  return out;
}

interface HomeIndex {
  byId: Map<string, Home>;
  byAlias: Map<string, Home>;
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
    errors.push({file: home.file, line: home.line, message: `duplicate id ~${home.id}, also the home of ${other.file}`});
    errors.push({file: other.file, line: other.line, message: `duplicate id ~${home.id}, also the home of ${home.file}`});
  }
  for (const home of homes)
    for (const alias of home.aliases) {
      const owner = byId.get(alias);
      if (owner) {
        errors.push({file: home.file, line: home.line, message: `alias ~${alias} is the id of ${owner.file}`});
        continue;
      }
      const other = byAlias.get(alias);
      if (other && other !== home) errors.push({file: home.file, line: home.line, message: `alias ~${alias} is also an alias of ${other.file}`});
      else byAlias.set(alias, home);
    }
  return {byId, byAlias};
}

class HomeChecker {
  private pathCache = new Map<string, string | null>();
  private stubIds = new Set<string>();

  constructor(private system: TypeSystem, private repoRoot: string, private index: HomeIndex, private set: HomeSet) {}

  check(home: Home, fm: Frontmatter): void {
    const error = (message: string, key?: string) => this.set.errors.push({file: home.file, line: key ? keyLine(fm, key) ?? home.line : home.line, message});
    const {type, data} = home;
    const isHelp = home.file.startsWith('public/help/');
    if (!home.name && (data.name === undefined || data.name === null))
      error('no name: set name:, title:, or start the body with a # heading');
    for (const [key, value] of Object.entries(data)) {
      if (key === 'feature' || key === 'id' || key === 'type' || value === null) continue;
      if (key === 'title' && typeof value === 'string') continue;
      if (key === 'part_of') {
        error('PART_OF is derived from the id path and never authored; remove part_of:', key);
        continue;
      }
      if (key === 'edges') {
        this.checkEdgesKey(home, value, error);
        continue;
      }
      const edge = this.system.keys.get(key);
      if (edge) {
        this.checkEdgeKey(home, key, edge, value, error);
        continue;
      }
      const member = type.members[key];
      if (member) {
        const problem = this.valueProblem(member, value);
        if (problem) error(`${key}: ${problem}`, key);
        continue;
      }
      if (isHelp && DOCUSAURUS_KEYS.includes(key)) continue;
      error(`unknown key '${key}' for type ${type.name}`, key);
    }
    for (const member of Object.values(type.members)) {
      if (member.nullable || member.name === 'id' || member.name === 'name') continue;
      if (data[member.name] === undefined || data[member.name] === null) error(`missing required key '${member.name}' for type ${type.name}`);
    }
    this.checkHierarchy(home, error);
    this.checkCitations(home, fm);
  }

  private checkHierarchy(home: Home, error: (message: string, key?: string) => void): void {
    const {type, localId} = home;
    if (!type.hierarchical) return;
    const segments = localId.split('/');
    if (type.root === 'feature' && segments.length === 2 && (home.data.owner === undefined || home.data.owner === null))
      error(`~${home.id} is a second-level feature node and must set owner:`);
    for (let i = 1; i < segments.length; i++) {
      const ancestor = (home.prefix ? `${home.prefix}:` : '') + segments.slice(0, i).join('/');
      if (this.index.byId.has(ancestor) || this.stubIds.has(ancestor)) continue;
      this.stubIds.add(ancestor);
      this.set.stubs.push({id: ancestor, type, neededBy: home.id});
      this.set.warnings.push({file: home.file, line: home.line, message: `no home for ~${ancestor} (parent of ~${home.id}); grok kg gen lists it as a stub`});
    }
  }

  private checkEdgeKey(home: Home, key: string, edge: EdgeType, value: unknown, error: (message: string, key?: string) => void): void {
    if (!Array.isArray(value)) {
      error(`${key}: expected a list of ${edge.name} targets, got ${JSON.stringify(value)}`, key);
      return;
    }
    const homeSide = edge.keySide;
    const otherSide = homeSide === 'from' ? 'to' : 'from';
    if (!edge[homeSide].some((t) => isSubtype(this.system, home.type.name, t))) {
      error(`${key}: ${edge.name}.${homeSide} must be ${edge[homeSide].map(pascal).join(' | ')}; this home is a ${home.type.name}`, key);
      return;
    }
    const targetKey = key === 'code' ? 'path' : 'to';
    value.forEach((item, i) => {
      const where = `${key}[${i}]`;
      let target: unknown;
      let props: Record<string, unknown> = {};
      if (typeof item === 'string') target = item;
      else if (item && typeof item === 'object' && !Array.isArray(item)) {
        ({[targetKey]: target, ...props} = item as Record<string, unknown>);
        if (target === undefined) {
          error(`${where}: a map item needs ${targetKey}:`, key);
          return;
        }
      } else {
        error(`${where}: expected an id or a map with ${targetKey}:, got ${JSON.stringify(item)}`, key);
        return;
      }
      if (typeof target !== 'string' || !target.trim()) {
        error(`${where}: ${targetKey} must be a string, got ${JSON.stringify(target)}`, key);
        return;
      }
      const problem = key === 'code' ? this.pathProblem(target) : this.refProblem(target, edge[otherSide]);
      if (problem) error(`${where}: ${problem}`, key);
      this.checkEdgeProperties(edge, props, where, error,key);
    });
  }

  private checkEdgesKey(home: Home, value: unknown, error: (message: string, key?: string) => void): void {
    if (!Array.isArray(value)) {
      error(`edges: expected a list of {type, to, ...} maps, got ${JSON.stringify(value)}`, 'edges');
      return;
    }
    value.forEach((item, i) => {
      const where = `edges[${i}]`;
      if (!item || typeof item !== 'object' || Array.isArray(item)) {
        error(`${where}: expected a map with type: and to:`, 'edges');
        return;
      }
      const {type, to, ...props} = item as Record<string, unknown>;
      const edge = typeof type === 'string' ? this.system.edges.get(type) : undefined;
      if (!edge) {
        error(`${where}: unknown edge type ${JSON.stringify(type)}`, 'edges');
        return;
      }
      if (edge.abstract || !edge.derivedBy.includes('annotation')) {
        error(`${where}: ${edge.name} is ${edge.abstract ? 'abstract' : 'never authored (derived_by lacks annotation)'}`, 'edges');
        return;
      }
      if (!edge.from.some((t) => isSubtype(this.system, home.type.name, t))) {
        error(`${where}: ${edge.name}.from must be ${edge.from.map(pascal).join(' | ')}; this home is a ${home.type.name}`, 'edges');
        return;
      }
      if (typeof to !== 'string' || !to.trim()) {
        error(`${where}: to must be an id, got ${JSON.stringify(to)}`, 'edges');
        return;
      }
      const problem = this.refProblem(to, edge.to);
      if (problem) error(`${where}: ${problem}`, 'edges');
      this.checkEdgeProperties(edge, props, where, error,'edges');
    });
  }

  private checkEdgeProperties(edge: EdgeType, props: Record<string, unknown>, where: string,
    error: (message: string, key?: string) => void, key: string): void {
    for (const [name, value] of Object.entries(props)) {
      const member = edge.properties[name];
      if (!member) {
        error(`${where}: ${edge.name} has no property '${name}' (${Object.keys(edge.properties).join(', ') || 'none'})`, key);
        continue;
      }
      if (value === null) continue;
      const problem = this.valueProblem(member, value);
      if (problem) error(`${where}.${name}: ${problem}`, key);
    }
    for (const member of Object.values(edge.properties))
      if (!member.nullable && (props[member.name] === undefined || props[member.name] === null))
        error(`${where}: ${edge.name} requires '${member.name}'`, key);
  }

  private valueProblem(member: Member, value: unknown): string | null {
    return checkValue(member, value, {
      provenance: this.system.provenance,
      path: (v) => this.pathProblem(v),
      ref: (v, m) => this.refProblem(v, m.refs!),
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
      problem = globSync(p, {cwd: this.repoRoot, posix: true, windowsPathsNoEscape: true}).length ? null : `no file matches '${value}'`;
    else
      problem = fs.existsSync(path.join(this.repoRoot, p)) ? null : `path '${value}' does not exist`;
    this.pathCache.set(value, problem);
    return problem;
  }

  /** Resolves [value] against the homes as one of [types] (kebab); counts extracted ids as unresolved-external. */
  refProblem(value: string, types: string[]): string | null {
    const raw = value.trim().replace(/^~/, '');
    const expanded = types.includes('node') ? [...this.system.nodes.keys()] : types;
    const prefixed = PREFIXED_ID.exec(raw);
    const candidates = prefixed ? [raw] : [...new Set(expanded.map((t) => {
      const prefix = this.system.nodes.get(t)?.prefix;
      return prefix ? `${prefix}:${raw}` : raw;
    }))];
    const hits = [...new Set(candidates.map((c) => this.index.byId.get(c) ?? this.index.byAlias.get(c)).filter((h): h is Home => !!h))];
    if (hits.length > 1) return `'${value}' is ambiguous: ${hits.map((h) => `~${h.id}`).join(', ')}; write the prefix`;
    if (hits.length === 1) {
      const home = hits[0];
      return types.some((t) => isSubtype(this.system, home.type.name, t)) ? null :
        `'${value}' resolves to ~${home.id}, a ${home.type.name}; expected ${types.map(pascal).join(' | ')}`;
    }
    const prefixType = prefixed ? this.system.prefixes.get(prefixed[1]) : undefined;
    const authoredTarget = prefixed ? (prefixType ? this.system.nodes.get(prefixType)!.authored : false) :
      expanded.some((t) => this.system.nodes.get(t)?.authored);
    if (!authoredTarget || (!prefixed && EXTERNAL_ID.test(raw))) {
      this.set.unresolvedExternal++;
      return null;
    }
    return `'${value}' does not resolve to any home document (as ${candidates.map((c) => `~${c}`).join(' or ')})`;
  }

  private checkCitations(home: Home, fm: Frontmatter): void {
    for (const {text, offset} of proseLines(home.body)) {
      const tokens: string[] = [];
      for (const m of text.matchAll(/`([^`\n]+)`/g)) tokens.push(m[1]);
      for (const m of text.matchAll(/\]\(([^)\s]+)\)/g)) tokens.push(m[1]);
      for (const token of tokens) {
        const cited = token.trim().replace(/#.*$/, '').replace(/:\d+(-\d+)?$/, '');
        if (!CITED_PATH.test(cited)) continue;
        const problem = this.pathProblem(cited);
        if (problem) this.set.errors.push({file: home.file, line: fm.bodyLine + offset, message: `cited ${problem}`});
      }
    }
  }
}

export function makeReport(system: TypeSystem, homes: HomeSet | null): Report {
  const byType: Record<string, number> = {};
  for (const home of homes?.homes ?? [])
    byType[home.type.name] = (byType[home.type.name] ?? 0) + 1;
  return {
    errors: [...system.errors, ...(homes?.errors ?? [])],
    warnings: [...system.warnings, ...(homes?.warnings ?? [])],
    types: {nodes: system.nodes.size, edges: system.edges.size, prefixes: system.prefixes.size},
    scanned: homes?.scanned ?? 0,
    homes: byType,
    unresolvedExternal: homes?.unresolvedExternal ?? 0,
    stubs: (homes?.stubs ?? []).map((s) => s.id),
  };
}
