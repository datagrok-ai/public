/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)

export type BytesKind = 'tables' | 'files';

export interface Ref {type: string; id?: string; nqName?: string}

export interface TypeSpec {
  route: string;
  saveRoute?: string;
  rank: number;
  tags?: boolean;
  listVia?: 'entities';
  typeId?: string;
  /** `afterSave` when the save may return another id and the bytes belong to that one. */
  bytes?: {kind: BytesKind; afterSave?: boolean; get(id: string): string; put(id: string): string};
  strip?(json: any): void;
  deps?(json: any): Ref[];
}

/** Extra listing rules a `--type` alias adds on top of its entity type. */
export interface TypeOptions {clause?: string; params?: Record<string, string>}

export const UUID = '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}';
export const UUID_RE = new RegExp(`^${UUID}$`, 'i');

export const isUuid = (token: string): boolean => UUID_RE.test(token);

export const nqNameOf = (json: any): string => `${json?.namespace ?? ''}${json?.name ?? json?.id ?? ''}`;

/** `/entities` records of some types (DataJob) carry no namespace; the server-side filter applied it. */
export const inNamespace = (json: any, namespace: string): boolean =>
  json?.namespace === undefined || json.namespace === null || json.namespace === namespace;

export function ref(type: string, id?: string): Ref {
  return {type, id};
}

export function idOnly(v: any): any {
  return v?.id ? {id: v.id} : v;
}

// Credentials.PASSWORD_PLACEHOLDER / SAME_PASSWORD (grok_shared/lib/src/credentials.dart)
export const PASSWORD_PLACEHOLDER = '_____________';
export const SAME_PASSWORD = 'not—changed';

/**
 * Never travels, whichever side asks: the platform's own connections, a space's `Files`
 * connection, and the personal `Home` share belong to the instance, not to the content.
 */
export function untransferableReason(type: string, json: any): {action: 'warn' | 'info'; reason: string} | null {
  if (type !== 'DataConnection') return null;
  if (json?.namespace === 'System:') return {action: 'info', reason: 'platform_connection'};
  if (json?.parameters?.isProject === true) return {action: 'info', reason: 'space_files_connection'};
  if (json?.name === 'Home') return {action: 'warn', reason: 'personal_storage'};
  return null;
}

export function stripCredentials(c: any): void {
  delete c.credentials;
  const params = c.parameters;
  if (!params) return;
  const masked = Object.keys(params).filter((k) => params[k] === PASSWORD_PLACEHOLDER || params[k] === SAME_PASSWORD);
  for (const k of masked)
    delete params[k];
  if (masked.length)
    c._credentials = masked;
}

const FILE_CALL_RE = /Open(?:ServerFile|File|Folder)[A-Za-z]*\s*\(\s*["']([^"']+)["']/g;

export function datasyncConnectionRefs(t: any): Ref[] {
  if (t.metaParams?.['.data-sync'] !== 'sync') return [];
  const script: string = t.metaParams?.['.script'] ?? '';
  const refs: Ref[] = [];
  for (const m of script.matchAll(FILE_CALL_RE)) {
    const nqName = m[1].split('/')[0];
    if (nqName)
      refs.push({type: 'DataConnection', nqName});
  }
  return refs;
}

function compactRelations(p: any): void {
  if (!Array.isArray(p.relations)) return;
  p.relations = p.relations
    // The space's own Files connection is stamped in by whichever server holds the space,
    // so its relation belongs to that server, not to the content.
    .filter((r: any) => r?.entity?.id && r.entity.parameters?.isProject !== true)
    .map((r: any) => ({id: r.id, entity: {'#type': 'EntityRecord', id: r.entity.id}, isLink: r.isLink ?? false}))
    .sort((a: any, b: any) => a.entity.id.localeCompare(b.entity.id));
}

export const TYPES: Record<string, TypeSpec> = {
  // POST is routed as `/groups/` (`routers/groups.dart:30`); `/groups` is a 404.
  UserGroup: {route: '/groups', saveRoute: '/groups/', rank: 0, strip: (g) => { g.parents = []; g.children = []; }},
  DataConnection: {route: '/connectors/connections', rank: 1, tags: true, strip: stripCredentials},
  Project: {route: '/projects', rank: 2, tags: true, strip: (p) => { delete p.storage; compactRelations(p); }},
  DataQuery: {
    route: '/connectors/queries', rank: 3, tags: true,
    strip: (q) => { q.connection = idOnly(q.connection); },
    deps: (q) => [ref('DataConnection', q.connection?.id)],
  },
  Script: {route: '/scripts', rank: 3, tags: true},
  TableInfo: {
    route: '/tables', rank: 4, tags: true,
    bytes: {kind: 'tables', get: (id) => `/tables/${id}/data`, put: (id) => `/tables/data?id=${id}`},
    strip: (t) => { delete t.fileInfo; delete t.columns; },
    deps: datasyncConnectionRefs,
  },
  ViewInfo: {
    route: '/views', rank: 5,
    strip: (v) => { v.table = idOnly(v.table); },
    deps: (v) => [ref('TableInfo', v.table?.id)],
  },
  ViewLayout: {route: '/layouts', rank: 5},
  FileInfo: {
    route: '/files', rank: 6, listVia: 'entities', typeId: '34d75630-e870-11e6-bfe1-590ff6f10d14',
    // `POST /files` dedups by (connection, dir, name) and answers with the existing row's
    // id (`files_service.dart:22-32`), so the bytes can only be written once it is known.
    bytes: {kind: 'files', afterSave: true, get: (id) => `/files/data/${id}`, put: (id) => `/files/data/${id}`},
    strip: (f) => { f.connection = idOnly(f.connection); delete f.tables; },
    deps: (f) => [ref('DataConnection', f.connection?.id)],
  },
  DataJob: {route: '/connectors/jobs', rank: 7, tags: true},
  Notebook: {
    route: '/notebooks', rank: 7, tags: true,
    // `environment` is a view of the ipynb kernelspec, and pushing it back overwrites the
    // kernel's display name with its id (`notebook.dart` `set environment`).
    strip: (n) => { n.tables = (n.tables ?? []).map(idOnly); delete n.environment; },
    deps: (n) => (n.tables ?? []).map((t: any) => ref('TableInfo', t?.id)),
  },
  PredictiveModelInfo: {
    route: '/ml', rank: 7, tags: true,
    strip: (m) => { m.trainedOn = idOnly(m.trainedOn); },
    deps: (m) => [ref('TableInfo', m.trainedOn?.id)],
  },
};

/** A root space is skipped by the default `/projects` listing — it has no namespace of its own. */
export const TYPE_OPTIONS: Record<string, TypeOptions> = {
  dashboard: {clause: 'isDashboard = true'},
  space: {clause: 'isDashboard = false and isEntity = false', params: {includeRoot: 'true'}},
};

export const TYPE_ALIASES: Record<string, string> = {
  conn: 'DataConnection', connection: 'DataConnection',
  query: 'DataQuery',
  script: 'Script',
  project: 'Project', dashboard: 'Project', space: 'Project',
  view: 'ViewInfo',
  layout: 'ViewLayout',
  table: 'TableInfo',
  file: 'FileInfo',
  group: 'UserGroup',
  job: 'DataJob',
  notebook: 'Notebook',
  model: 'PredictiveModelInfo',
};

export const DEFAULT_TYPES: string[] = Object.keys(TYPES);

export function typeOptionsOf(alias: string): TypeOptions | undefined {
  return TYPE_OPTIONS[alias.trim().toLowerCase()];
}

export function resolveType(alias: string): string {
  const name = TYPE_ALIASES[alias.trim().toLowerCase()] ?? alias.trim();
  if (!TYPES[name])
    throw new Error(`Entity type '${alias}' is not supported. Supported: ${Object.keys(TYPES).join(', ')}`);
  return name;
}

/** One entity type can be asked for under several aliases, but only under one listing rule. */
export function resolveTypes(aliases: string[]): {types: string[]; typeOptions: Record<string, TypeOptions>} {
  const types: string[] = [];
  const typeOptions: Record<string, TypeOptions> = {};
  const aliasOf: Record<string, string> = {};
  for (const alias of aliases) {
    const type = resolveType(alias);
    const options = typeOptionsOf(alias);
    if (types.includes(type)) {
      if (JSON.stringify(typeOptions[type] ?? null) !== JSON.stringify(options ?? null))
        throw new Error(`--type ${aliasOf[type]} and ${alias.trim()} both select ${type} but list it differently — pass only one`);
      continue;
    }
    types.push(type);
    aliasOf[type] = alias.trim();
    if (options)
      typeOptions[type] = options;
  }
  return {types, typeOptions};
}

export function rankOf(type: string): number {
  return TYPES[type]?.rank ?? 99;
}

/**
 * `Chem:Chembl` → `Chem.Chembl`; FileInfo `Admin:Data/a/b.csv` → `Admin.Data.a.b.csv`.
 * A namespace-less entity is named by its full path, so two files with the same basename
 * in different directories do not claim the same bundle file.
 */
export function fileNameFor(json: any): string {
  const name: string = json?.name ?? '';
  if (!name) return json?.id ?? 'unnamed';
  const nqName = json.namespace ? `${json.namespace}${name}` : (json.path ?? name);
  const cleaned = String(nqName).replace(/[:/\\?*"<>|]+/g, '.').replace(/^\.+|\.+$/g, '');
  // Windows resolves `NUL.json` (and `COM1.csv.json`) to a device, whatever the extension.
  const head = cleaned.split('.')[0];
  return RESERVED_RE.test(head) ? `${head}-${json.id ?? 'entity'}${cleaned.slice(head.length)}` : cleaned;
}

const RESERVED_RE = /^(NUL|CON|PRN|AUX|COM[1-9]|LPT[1-9])$/i;
