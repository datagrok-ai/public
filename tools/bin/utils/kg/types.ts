/// The knowledge-graph type system: `schema.yaml`, `nodes/*.yaml`, `edges/*.yaml`
/// (core/docs/knowledge-graph/CONVENTIONS.md §7). Loads the files, parses every member
/// with the TypeScript compiler and checks the constraints listed in `schema.yaml`.
import * as fs from 'fs';
import * as path from 'path';
import * as yaml from 'js-yaml';
import {ts} from 'ts-morph';

export type Scalar = 'string' | 'number' | 'boolean' | 'Date' | 'Text' | 'Version' | 'Path' | 'Url' |
  'Provenance' | 'TypeUnion' | 'Member';

export interface Member {
  name: string;
  nullable: boolean;
  list: boolean;
  kind: 'scalar' | 'enum' | 'ref' | 'record';
  spec: string;
  scalar?: Scalar;
  literals?: string[];
  /** Kebab node type names; `node` stands for `Node`, which admits anything. */
  refs?: string[];
  recordValue?: 'string' | 'Member';
  default?: string | number | boolean;
}

export interface Issue {
  file: string;
  line?: number;
  message: string;
}

export interface NodeType {
  name: string;
  file: string;
  extends?: string;
  /** Self first, `node` last. */
  chain: string[];
  root?: string;
  prefix?: string;
  abstract: boolean;
  hierarchical: boolean;
  authored: boolean;
  visibility?: string;
  was?: string;
  description: string;
  home?: string;
  own: Record<string, Member>;
  members: Record<string, Member>;
  inherit: string[];
}

export interface EdgeType {
  name: string;
  file: string;
  extends?: string;
  chain: string[];
  abstract: boolean;
  was?: string;
  from: string[];
  to: string[];
  key?: string;
  keySide: 'from' | 'to';
  cardinality: 'one' | 'many';
  sameType: boolean;
  symmetric: boolean;
  acyclic: boolean;
  derivedBy: string[];
  description: string;
  own: Record<string, Member>;
  properties: Record<string, Member>;
}

export interface TypeSystem {
  nodes: Map<string, NodeType>;
  edges: Map<string, EdgeType>;
  /** Prefix -> node type that declares it. */
  prefixes: Map<string, string>;
  /** Edge key -> edge type. */
  keys: Map<string, EdgeType>;
  roots: string[];
  provenance: string[];
  reservedNodeFields: string[];
  errors: Issue[];
  warnings: Issue[];
}

export interface ParseOptions {
  /** Kebab node type names a reference may name; `node` is always known. */
  types?: Iterable<string>;
  /** Inside schema.yaml `Record<string, Member>` is allowed as well. */
  schema?: boolean;
}

const SCALARS: Scalar[] = ['string', 'number', 'boolean', 'Date', 'Text', 'Version', 'Path', 'Url',
  'Provenance', 'TypeUnion', 'Member'];
const MEMBER_NAME = /^[a-z][a-z0-9_]*$/;
const NODE_TYPE_NAME = /^[a-z][a-z0-9-]*$/;
const EDGE_TYPE_NAME = /^[A-Z][A-Z0-9_]*$/;
const PREFIX = /^[A-Z][A-Za-z]{0,5}$/;
const DEFAULT_ROOTS = ['feature', 'concept', 'component', 'artifact', 'work', 'actor', 'infra', 'type'];
const DEFAULT_PROVENANCE = ['annotation', 'filesystem', 'ast', 'registry', 'git', 'external', 'llm', 'manual'];
const DEFAULT_RESERVED = ['id', 'name', 'description', 'status', 'visibility', 'owner', 'aliases',
  'source_layer', 'home', 'provenance', 'batch'];
const EDGE_BUILD_FIELDS = ['derived_by', 'confidence', 'evidence', 'batch'];
const MAX_DEPTH = 3;

export function pascal(kebab: string): string {
  return kebab.split(/[-_]/).map((w) => w.charAt(0).toUpperCase() + w.slice(1).toLowerCase()).join('');
}

export function firstSentence(text: string): string {
  const flat = String(text ?? '').split(/\s+/).join(' ').trim();
  const m = /(.+?\.)(\s|$)/.exec(flat);
  return m ? m[1] : flat;
}

class MemberError extends Error {}

/** Parses one `<name>[?]: <ts-type> [= <default>]` entry (CONVENTIONS §7.1). */
export function parseMember(name: string, spec: unknown, options: ParseOptions = {}): {member?: Member, error?: string} {
  try {
    return {member: parseMemberOrThrow(name, spec, options)};
  } catch (e: any) {
    if (e instanceof MemberError) return {error: `${name.replace(/\?$/, '')}: ${e.message}`};
    throw e;
  }
}

function parseMemberOrThrow(rawName: string, spec: unknown, options: ParseOptions): Member {
  if (typeof spec !== 'string')
    throw new MemberError(`the type must be a string, got ${JSON.stringify(spec)}`);
  const nullable = rawName.endsWith('?');
  const name = nullable ? rawName.slice(0, -1) : rawName;
  if (!MEMBER_NAME.test(name))
    throw new MemberError(`bad member name (expected ${MEMBER_NAME})`);
  const source = ts.createSourceFile('member.ts', `class _ { ${rawName}: ${spec}; }`, ts.ScriptTarget.Latest, true, ts.ScriptKind.TS);
  const diagnostics: ts.Diagnostic[] = (source as any).parseDiagnostics ?? [];
  if (diagnostics.length)
    throw new MemberError(`syntax error in '${spec}': ${ts.flattenDiagnosticMessageText(diagnostics[0].messageText, ' ')}`);
  const cls = source.statements[0];
  if (source.statements.length !== 1 || !ts.isClassDeclaration(cls) || cls.members.length !== 1)
    throw new MemberError(`'${spec}' is not a single type expression`);
  const prop = cls.members[0];
  if (!ts.isPropertyDeclaration(prop) || !prop.type)
    throw new MemberError(`'${spec}' is not a type`);
  const known = new Set<string>(['node', ...(options.types ?? [])]);
  const member = readType(prop.type, known, !!options.schema, false);
  member.name = name;
  member.nullable = nullable;
  member.spec = spec;
  if (prop.initializer) {
    if (!nullable) throw new MemberError('a default is allowed only on a nullable member');
    member.default = readDefault(prop.initializer);
    const problem = defaultProblem(member);
    if (problem) throw new MemberError(problem);
  }
  return member;
}

function readType(node: ts.TypeNode, known: Set<string>, schema: boolean, inList: boolean): Member {
  const base: Member = {name: '', nullable: false, list: false, kind: 'scalar', spec: ''};
  switch (node.kind) {
    case ts.SyntaxKind.StringKeyword: return {...base, scalar: 'string'};
    case ts.SyntaxKind.NumberKeyword: return {...base, scalar: 'number'};
    case ts.SyntaxKind.BooleanKeyword: return {...base, scalar: 'boolean'};
  }
  if (ts.isParenthesizedTypeNode(node))
    return readType(node.type, known, schema, inList);
  if (ts.isArrayTypeNode(node)) {
    if (inList) throw new MemberError('nested arrays are not allowed');
    const element = readType(node.elementType, known, schema, true);
    if (element.kind === 'record') throw new MemberError('an array of Record is not allowed');
    return {...element, list: true};
  }
  if (ts.isLiteralTypeNode(node)) {
    if (!ts.isStringLiteral(node.literal)) throw new MemberError(`only string literals may form an enum, got '${node.literal.getText()}'`);
    return {...base, kind: 'enum', literals: [node.literal.text]};
  }
  if (ts.isUnionTypeNode(node)) {
    const parts = node.types.map((t) => readType(t, known, schema, inList));
    if (parts.every((p) => p.kind === 'enum'))
      return {...base, kind: 'enum', literals: parts.flatMap((p) => p.literals!)};
    if (parts.every((p) => p.kind === 'ref'))
      return {...base, kind: 'ref', refs: [...new Set(parts.flatMap((p) => p.refs!))]};
    throw new MemberError(`a union must be all string literals or all node types, got '${node.getText()}'`);
  }
  if (ts.isTypeReferenceNode(node)) {
    if (!ts.isIdentifier(node.typeName)) throw new MemberError(`qualified type '${node.getText()}' is not allowed`);
    const typeName = node.typeName.text;
    if (typeName === 'Record') {
      const args = (node.typeArguments ?? []).map((a) => a.getText());
      if (args.join(', ') === 'string, string') return {...base, kind: 'record', recordValue: 'string'};
      if (schema && args.join(', ') === 'string, Member') return {...base, kind: 'record', recordValue: 'Member'};
      throw new MemberError(`only Record<string, string> is allowed, got '${node.getText()}'`);
    }
    if (node.typeArguments) throw new MemberError(`generic type '${node.getText()}' is not allowed`);
    if ((SCALARS as string[]).includes(typeName)) return {...base, scalar: typeName as Scalar};
    const kebab = typeName === 'Node' ? 'node' : [...known].find((t) => pascal(t) === typeName);
    if (!kebab || !known.has(kebab)) throw new MemberError(`unknown type '${typeName}'`);
    return {...base, kind: 'ref', refs: [kebab]};
  }
  throw new MemberError(`'${node.getText()}' (${ts.SyntaxKind[node.kind]}) is not allowed here`);
}

function readDefault(node: ts.Expression): string | number | boolean {
  if (ts.isStringLiteral(node)) return node.text;
  if (ts.isNumericLiteral(node)) return Number(node.text);
  if (node.kind === ts.SyntaxKind.TrueKeyword) return true;
  if (node.kind === ts.SyntaxKind.FalseKeyword) return false;
  if (ts.isPrefixUnaryExpression(node) && node.operator === ts.SyntaxKind.MinusToken && ts.isNumericLiteral(node.operand))
    return -Number(node.operand.text);
  throw new MemberError(`a default must be a string, number or boolean literal, got '${node.getText()}'`);
}

function defaultProblem(m: Member): string | null {
  const d = m.default;
  if (m.list || m.kind === 'ref' || m.kind === 'record') return `a default is not allowed on a ${m.list ? 'list' : m.kind}`;
  if (m.kind === 'enum') return m.literals!.includes(String(d)) && typeof d === 'string' ? null : `default '${d}' is not one of ${m.literals!.map((l) => `'${l}'`).join(' | ')}`;
  const expected = m.scalar === 'number' ? 'number' : m.scalar === 'boolean' ? 'boolean' : 'string';
  return typeof d === expected ? null : `default ${JSON.stringify(d)} is not a ${m.scalar}`;
}

/** Parses a `Feature | Ticket` value into kebab node type names. */
export function parseTypeUnion(spec: unknown, known: Iterable<string>): {types?: string[], error?: string} {
  if (typeof spec !== 'string' || !spec.trim()) return {error: `expected a union of node type names, got ${JSON.stringify(spec)}`};
  const byPascal = new Map<string, string>([['Node', 'node']]);
  for (const t of known) byPascal.set(pascal(t), t);
  const types: string[] = [];
  for (const part of spec.split('|').map((s) => s.trim())) {
    const kebab = byPascal.get(part);
    if (!kebab) return {error: `unknown node type '${part}' in '${spec}'`};
    if (!types.includes(kebab)) types.push(kebab);
  }
  return {types};
}

export function isSubtype(system: TypeSystem, type: string, ancestor: string): boolean {
  if (ancestor === 'node' || type === ancestor) return true;
  return system.nodes.get(type)?.chain.includes(ancestor) ?? false;
}

/** Why [child] is not a narrowing of [anc]; null when it is. */
export function narrowingProblem(child: Member, anc: Member, subtype: (t: string, a: string) => boolean): string | null {
  if (child.list !== anc.list) return `changes ${anc.list ? 'a list into a single value' : 'a single value into a list'}`;
  if (child.kind !== anc.kind) return `changes kind from ${describeKind(anc)} to ${describeKind(child)}`;
  if (!anc.nullable && child.nullable) return 'loosens a required member to nullable';
  switch (child.kind) {
    case 'scalar':
      return child.scalar === anc.scalar ? null : `changes ${anc.scalar} to ${child.scalar}`;
    case 'enum': {
      const extra = child.literals!.filter((l) => !anc.literals!.includes(l));
      return extra.length ? `adds ${extra.map((l) => `'${l}'`).join(', ')} not allowed by the ancestor` : null;
    }
    case 'ref': {
      if (anc.refs!.includes('node')) return null;
      const extra = child.refs!.filter((r) => !anc.refs!.some((a) => subtype(r, a)));
      return extra.length ? `references ${extra.map(pascal).join(' | ')}, not a subtype of ${anc.refs!.map(pascal).join(' | ')}` : null;
    }
    case 'record':
      return child.recordValue === anc.recordValue ? null : 'changes the Record value type';
  }
}

function describeKind(m: Member): string {
  return m.kind === 'scalar' ? m.scalar! : m.kind === 'ref' ? 'reference' : m.kind;
}

interface RawType {
  file: string;
  data: Record<string, any>;
}

export function loadTypeSystem(kgRoot: string): TypeSystem {
  const system: TypeSystem = {
    nodes: new Map(), edges: new Map(), prefixes: new Map(), keys: new Map(),
    roots: DEFAULT_ROOTS, provenance: DEFAULT_PROVENANCE, reservedNodeFields: DEFAULT_RESERVED,
    errors: [], warnings: [],
  };
  const error = (file: string, message: string, line?: number) => system.errors.push({file, line, message});
  const schemaFile = path.join(kgRoot, 'schema.yaml');
  const schema = readYaml(schemaFile, error);
  if (!schema) return system;
  const constraints = schema.constraints ?? {};
  if (Array.isArray(constraints.roots)) system.roots = constraints.roots.map(String);
  const reserved = schema.manifest?.reserved_fields;
  if (Array.isArray(reserved)) system.reservedNodeFields = reserved.map(String);
  const provenanceSpec = parseMember('provenance', schema.member?.aliases?.Provenance, {schema: true});
  if (provenanceSpec.member?.kind === 'enum') system.provenance = provenanceSpec.member.literals!;
  else error(schemaFile, `member.aliases.Provenance: ${provenanceSpec.error ?? 'expected a string literal union'}`);

  const rawNodes = loadFolder(path.join(kgRoot, 'nodes'), NODE_TYPE_NAME, 'node', error);
  const rawEdges = loadFolder(path.join(kgRoot, 'edges'), EDGE_TYPE_NAME, 'edge', error);
  // a file that failed to load still names its type (nodes/<type>.yaml), so references to it do not cascade into noise
  const nodeNames = [...rawNodes.keys(), ...unreadableTypes(path.join(kgRoot, 'nodes'), rawNodes)];
  const schemaNode = parseSchemaSection(schema.node_type, `${schemaFile} node_type`, error);
  const schemaEdge = parseSchemaSection(schema.edge_type, `${schemaFile} edge_type`, error);
  for (const raw of rawNodes.values()) checkTopLevel(raw, schemaNode, system, error);
  for (const raw of rawEdges.values()) checkTopLevel(raw, schemaEdge, system, error);

  const nodeChains = buildChains(rawNodes, error);
  const edgeChains = buildChains(rawEdges, error);
  for (const [name, raw] of rawNodes) {
    const chain = nodeChains.get(name)!;
    const d = raw.data;
    const inherited = (key: string) => chain.map((t) => rawNodes.get(t)?.data[key]).find((v) => v !== undefined);
    const own = parseProperties(raw, nodeNames, error);
    const isBase = name === (constraints.base ?? 'node');
    const last = rawNodes.get(chain[chain.length - 1]);
    if (last && last.data.extends === undefined && chain[chain.length - 1] !== 'node')
      error(raw.file, `node ${name}: the extends chain does not end at node: ${chain.join(' -> ')}`);
    if (!isBase && d.extends === undefined) error(raw.file, `node ${name}: no extends`);
    if (system.roots.includes(name) && d.extends !== 'node') error(raw.file, `root ${name} must extend node directly`);
    if (!isBase && d.extends === 'node' && !system.roots.includes(name))
      error(raw.file, `node ${name}: extends node directly but is not one of the roots (${system.roots.join(', ')})`);
    if (chain.length - 2 > MAX_DEPTH) error(raw.file, `node ${name}: ${chain.length - 2} levels below the root, at most ${MAX_DEPTH} allowed (${chain.join(' -> ')})`);
    if (!isBase)
      for (const m of Object.keys(own))
        if (system.reservedNodeFields.includes(m)) error(raw.file, `node ${name}.${m}: reserved field name, declared by the base or written by the build`);
    system.nodes.set(name, {
      name, file: raw.file, extends: d.extends, chain,
      root: chain.length >= 2 ? chain[chain.length - 2] : undefined,
      prefix: inherited('prefix'),
      abstract: d.abstract === true,
      hierarchical: inherited('hierarchical') === true,
      authored: inherited('authored') === true,
      visibility: inherited('visibility'),
      was: d.was, description: String(d.description ?? ''), home: d.home === undefined ? undefined : String(d.home),
      own, members: {}, inherit: Array.isArray(d.inherit) ? d.inherit.map(String) : [],
    });
  }
  for (const missing of system.roots.filter((r) => !rawNodes.has(r)))
    error(path.join(kgRoot, 'nodes'), `root ${missing} is missing`);
  for (const node of system.nodes.values())
    node.members = mergeMembers(node.chain.map((t) => system.nodes.get(t)?.own ?? {}));
  checkPrefixes(system, rawNodes, error);
  const subtype = (t: string, a: string) => isSubtype(system, t, a);
  for (const node of system.nodes.values())
    checkMemberNarrowing(node.name, node.own, node.chain.slice(1).map((t) => system.nodes.get(t)), node.file, subtype, error);

  for (const [name, raw] of rawEdges) {
    const d = raw.data;
    const chain = edgeChains.get(name)!;
    const from = parseTypeUnion(d.from, nodeNames);
    const to = parseTypeUnion(d.to, nodeNames);
    if (from.error) error(raw.file, `edge ${name}.from: ${from.error}`);
    if (to.error) error(raw.file, `edge ${name}.to: ${to.error}`);
    const own = parseProperties(raw, nodeNames, error);
    for (const m of Object.keys(own))
      if (EDGE_BUILD_FIELDS.includes(m)) error(raw.file, `edge ${name}.${m}: reserved field name, written by the build`);
    const derivedBy = Array.isArray(d.derived_by) ? d.derived_by.map(String) : [];
    if (d.key !== undefined && !derivedBy.includes('annotation'))
      error(raw.file, `edge ${name}: has key '${d.key}' but 'annotation' is not in derived_by`);
    if (d.same_type === true && from.types && to.types && [...from.types].sort().join('|') !== [...to.types].sort().join('|'))
      error(raw.file, `edge ${name}: same_type but from (${d.from}) differs from to (${d.to})`);
    system.edges.set(name, {
      name, file: raw.file, extends: d.extends, chain, abstract: d.abstract === true, was: d.was,
      from: from.types ?? [], to: to.types ?? [],
      key: d.key === undefined ? undefined : String(d.key),
      keySide: d.key_side === 'to' ? 'to' : 'from',
      cardinality: d.cardinality === 'one' ? 'one' : 'many',
      sameType: d.same_type === true, symmetric: d.symmetric === true, acyclic: d.acyclic === true,
      derivedBy, description: String(d.description ?? ''), own, properties: {},
    });
  }
  for (const edge of system.edges.values()) {
    edge.properties = mergeMembers(edge.chain.map((t) => system.edges.get(t)?.own ?? {}));
    if (edge.key !== undefined) {
      const other = system.keys.get(edge.key);
      if (other) error(edge.file, `edge key '${edge.key}' is used by both ${other.name} and ${edge.name}`);
      else system.keys.set(edge.key, edge);
    }
    const ancestors = edge.chain.slice(1).map((t) => system.edges.get(t));
    checkMemberNarrowing(edge.name, edge.own, ancestors, edge.file, subtype, error);
    for (const side of ['from', 'to'] as const)
      for (const anc of ancestors) {
        if (!anc || anc[side].includes('node')) continue;
        const outside = edge[side].filter((t) => !anc[side].some((a) => subtype(t, a)));
        if (!outside.length) continue;
        error(edge.file, `edge ${edge.name}.${side}: ${outside.map(pascal).join(' | ')} is not within ${anc.name}.${side}: ${anc[side].map(pascal).join(' | ')}`);
        break;
      }
  }
  return system;
}

function readYaml(file: string, error: (file: string, message: string, line?: number) => void): Record<string, any> | null {
  let text: string;
  try {
    text = fs.readFileSync(file, 'utf8');
  } catch (e: any) {
    error(file, `cannot read: ${e.message}`);
    return null;
  }
  try {
    const data = yaml.load(text, {filename: file});
    if (data === null || typeof data !== 'object' || Array.isArray(data)) {
      error(file, 'expected a YAML mapping');
      return null;
    }
    return data as Record<string, any>;
  } catch (e: any) {
    error(file, `YAML error: ${e.reason ?? e.message}`, e.mark ? e.mark.line + 1 : undefined);
    return null;
  }
}

function loadFolder(folder: string, namePattern: RegExp, kind: string,
  error: (file: string, message: string, line?: number) => void): Map<string, RawType> {
  const out = new Map<string, RawType>();
  if (!fs.existsSync(folder)) {
    error(folder, 'folder not found');
    return out;
  }
  for (const entry of fs.readdirSync(folder).filter((f) => f.endsWith('.yaml')).sort()) {
    const file = path.join(folder, entry);
    const data = readYaml(file, error);
    if (!data) continue;
    if (typeof data.type !== 'string') {
      error(file, 'no type key');
      continue;
    }
    const name: string = data.type;
    if (!namePattern.test(name)) error(file, `${kind} type name '${name}' does not match ${namePattern}`);
    if (entry !== `${name}.yaml`) error(file, `file name does not match type '${name}' (expected ${name}.yaml)`);
    if (out.has(name)) {
      error(file, `duplicate ${kind} type '${name}', also declared in ${out.get(name)!.file}`);
      continue;
    }
    out.set(name, {file, data});
  }
  return out;
}

function unreadableTypes(folder: string, loaded: Map<string, RawType>): string[] {
  if (!fs.existsSync(folder)) return [];
  const loadedFiles = new Set([...loaded.values()].map((r) => path.basename(r.file)));
  return fs.readdirSync(folder).filter((f) => f.endsWith('.yaml') && !loadedFiles.has(f)).map((f) => f.slice(0, -5));
}

function parseSchemaSection(section: unknown, where: string,
  error: (file: string, message: string) => void): Record<string, Member> {
  const members: Record<string, Member> = {};
  if (!section || typeof section !== 'object') {
    error(where, 'schema section is missing');
    return members;
  }
  for (const [name, spec] of Object.entries(section as Record<string, unknown>)) {
    const parsed = parseMember(name, spec, {schema: true});
    if (parsed.member) members[parsed.member.name] = parsed.member;
    else error(where, parsed.error!);
  }
  return members;
}

/** Checks a type file's top-level keys against the schema's `node_type` / `edge_type` members. */
function checkTopLevel(raw: RawType, schema: Record<string, Member>, system: TypeSystem,
  error: (file: string, message: string) => void): void {
  const d = raw.data;
  for (const key of Object.keys(d))
    if (!schema[key]) error(raw.file, `${d.type}: unknown key '${key}' (allowed: ${Object.keys(schema).join(', ')})`);
  for (const m of Object.values(schema)) {
    const value = d[m.name];
    if (value === undefined || value === null) {
      if (!m.nullable) error(raw.file, `${d.type}: missing required key '${m.name}'`);
      continue;
    }
    const problem = checkValue(m, value, {provenance: system.provenance});
    if (problem) error(raw.file, `${d.type}.${m.name}: ${problem}`);
  }
}

export interface ValueHooks {
  provenance: string[];
  /** Returns a problem for a Path value, or null. */
  path?: (value: string) => string | null;
  /** Returns a problem for a reference value, or null. */
  ref?: (value: string, member: Member) => string | null;
}

/** Why [value] does not fit [member]; null when it does. Reference and Path checks go through [hooks]. */
export function checkValue(member: Member, value: unknown, hooks: ValueHooks): string | null {
  if (member.list) {
    if (!Array.isArray(value)) return `expected a list, got ${JSON.stringify(value)}`;
    for (const item of value) {
      const problem = checkValue({...member, list: false}, item, hooks);
      if (problem) return problem;
    }
    return null;
  }
  const shown = JSON.stringify(value instanceof Date ? value.toISOString() : value);
  switch (member.kind) {
    case 'enum':
      return typeof value === 'string' && member.literals!.includes(value) ? null :
        `${shown} is not one of ${member.literals!.map((l) => `'${l}'`).join(' | ')}`;
    case 'ref':
      if (typeof value !== 'string' || !value.trim()) return `expected an id, got ${shown}`;
      return hooks.ref ? hooks.ref(value, member) : null;
    case 'record':
      if (typeof value !== 'object' || value === null || Array.isArray(value)) return `expected a mapping, got ${shown}`;
      for (const [k, v] of Object.entries(value))
        if (typeof v !== 'string') return `${k}: expected a string, got ${JSON.stringify(v)}`;
      return null;
  }
  switch (member.scalar) {
    case 'number': return typeof value === 'number' ? null : `expected a number, got ${shown}`;
    case 'boolean': return typeof value === 'boolean' ? null : `expected a boolean, got ${shown}`;
    case 'Date':
      if (value instanceof Date) return null;
      return typeof value === 'string' && /^\d{4}-\d{2}(-\d{2}([T ]\d{2}:\d{2}(:\d{2}(\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?)?)?$/.test(value) ? null :
        `expected an ISO date or datetime, got ${shown}`;
    case 'Version':
      if (typeof value === 'number') return 'a Version must be a quoted string; YAML read the unquoted value as a number, quote the value';
      return typeof value === 'string' && /^\d+(\.\d+)*$/.test(value) ? null : `expected a version like 1.28.0, got ${shown}`;
    case 'Url':
      return typeof value === 'string' && /^[a-z][a-z0-9+.-]*:\/\//i.test(value) ? null : `expected an absolute URL, got ${shown}`;
    case 'Provenance':
      return typeof value === 'string' && hooks.provenance.includes(value) ? null :
        `${shown} is not one of ${hooks.provenance.map((l) => `'${l}'`).join(' | ')}`;
    case 'Path':
      if (typeof value !== 'string' || !value.trim()) return `expected a path, got ${shown}`;
      return hooks.path ? hooks.path(value) : null;
    default:
      return typeof value === 'string' ? null : `expected a string, got ${shown}`;
  }
}

function buildChains(raw: Map<string, RawType>, error: (file: string, message: string) => void): Map<string, string[]> {
  const chains = new Map<string, string[]>();
  for (const [name, r] of raw) {
    const chain: string[] = [];
    let t: string | undefined = name;
    while (t !== undefined) {
      if (chain.includes(t)) {
        error(r.file, `cycle in extends at ${t}: ${[...chain, t].join(' -> ')}`);
        break;
      }
      chain.push(t);
      const next: RawType | undefined = raw.get(t);
      if (!next) {
        error(r.file, `${name}: extends target '${t}' does not exist`);
        break;
      }
      t = next.data.extends === undefined ? undefined : String(next.data.extends);
    }
    chains.set(name, chain);
  }
  return chains;
}

function parseProperties(raw: RawType, nodeNames: string[], error: (file: string, message: string) => void): Record<string, Member> {
  const out: Record<string, Member> = {};
  const props = raw.data.properties;
  if (props === undefined || props === null) return out;
  if (typeof props !== 'object' || Array.isArray(props)) {
    error(raw.file, `${raw.data.type}: properties must be a mapping`);
    return out;
  }
  for (const [name, spec] of Object.entries(props)) {
    const parsed = parseMember(name, spec, {types: nodeNames});
    if (parsed.member) out[parsed.member.name] = parsed.member;
    else error(raw.file, `${raw.data.type}.${parsed.error}`);
  }
  return out;
}

/** Ancestors first, so a child's declaration replaces its ancestor's. */
function mergeMembers(chainOwn: Record<string, Member>[]): Record<string, Member> {
  const out: Record<string, Member> = {};
  for (const own of [...chainOwn].reverse())
    Object.assign(out, own);
  return out;
}

function checkPrefixes(system: TypeSystem, raw: Map<string, RawType>, error: (file: string, message: string) => void): void {
  for (const [name, r] of raw) {
    const prefix = r.data.prefix;
    if (prefix === undefined) continue;
    if (!PREFIX.test(String(prefix))) error(r.file, `node ${name}: bad prefix '${prefix}' (expected ${PREFIX})`);
    const owner = system.prefixes.get(String(prefix));
    if (owner) error(r.file, `prefix ${prefix} is declared by both ${owner} and ${name}`);
    else system.prefixes.set(String(prefix), name);
    for (const anc of system.nodes.get(name)!.chain.slice(1))
      if (raw.get(anc)?.data.prefix !== undefined) error(r.file, `node ${name}: redeclares prefix inherited from ${anc}`);
  }
}

function checkMemberNarrowing(typeName: string, own: Record<string, Member>, ancestors: ({own: Record<string, Member>, name: string} | undefined)[],
  file: string, subtype: (t: string, a: string) => boolean, error: (file: string, message: string) => void): void {
  for (const [name, member] of Object.entries(own))
    for (const anc of ancestors) {
      const theirs = anc?.own[name];
      if (!theirs) continue;
      const problem = narrowingProblem(member, theirs, subtype);
      if (!problem) continue;
      error(file, `${typeName}.${name}: ${problem} (declared by ${anc!.name} as ${theirs.spec})`);
      break;
    }
}

/** Node types in glossary order: the base, then each root's subtree by depth and name. */
export function nodeOrder(system: TypeSystem): NodeType[] {
  const order = ['node', ...system.roots];
  const key = (n: NodeType) => {
    const index = order.indexOf(n.name === 'node' ? 'node' : n.root ?? '');
    return [index < 0 ? 99 : index, n.chain.length, n.name] as const;
  };
  return [...system.nodes.values()].sort((a, b) => {
    const ka = key(a), kb = key(b);
    return ka[0] - kb[0] || ka[1] - kb[1] || (ka[2] < kb[2] ? -1 : ka[2] > kb[2] ? 1 : 0);
  });
}

/** Edge types in glossary order: grouped by the base they end at, then by depth and name. */
export function edgeOrder(system: TypeSystem): EdgeType[] {
  return [...system.edges.values()].sort((a, b) => {
    const ba = a.chain[a.chain.length - 1], bb = b.chain[b.chain.length - 1];
    if (ba !== bb) return ba < bb ? -1 : 1;
    return a.chain.length - b.chain.length || (a.name < b.name ? -1 : a.name > b.name ? 1 : 0);
  });
}
