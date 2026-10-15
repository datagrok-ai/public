/**
 * Typed surface for entity-mapped domain tables (`grok.dapi.domains`): condition-tree
 * filter types, typed results, the DomainError class family, predicate helpers, the fluent
 * query builder, and optimistic-concurrency helpers.
 *
 * Two rules hold everywhere: condition VALUES are bound server-side — never interpolated
 * into filter strings, so any string value is expressible (apostrophes included, which the
 * smart-filter string grammar cannot quote); and datetime columns materialize as dayjs on
 * JSON reads — every client, typed or not, resolves them from the domain registry
 * (`DomainTableClientOptions` overrides exist for callers that must avoid it).
 *
 * @remarks BREAKING (codegen v2, GROK-20602): `grok api`-generated clients type datetime
 * columns as `Dayjs` and thread `<Table>Column`/`<Table>Expand` generics; regenerating
 * db.ts changes its surface. The generated `dapi2.domains` namespace is removed
 * (GROK-20601) — this module and `DomainTableClient` are the API.
 */
import type {Dayjs} from 'dayjs';
import type {DataFrame} from './dataframe';

// The runtime DG namespace, declared instead of imported: dg.ts re-exports this module,
// so importing it back would close a module cycle (same seam as IDomainQueryExecutor).
declare let DG: any;

/** The system columns every domain table carries (always projected on reads). */
export type DomainSystemColumn = 'id' | 'version' | 'created_on' | 'updated_on' | 'author_id';

/** Column reference: a declared column of the table, a system column, or a dotted FK path
 * ('project_id.name', up to 3 forward hops). Used by the condition, sort, and aggregate
 * column slots — those accept system columns even when a hand-written [TColumn] union omits
 * them (generated `<Table>Column` unions include them). Other slots (equality-map `where`,
 * `columns`, `select()`, `fetchFields()`, column security) stay TColumn-typed. */
export type DomainColumnRef<TColumn extends string = string> =
  TColumn | DomainSystemColumn | `${string}.${string}`;

/** Operators of the canonical condition tree (server: filter_compiler.dart:149-153 + fuzzy). */
export type DomainConditionOperator =
  '=' | '!=' | '>' | '>=' | '<' | '<=' | 'like' | 'not like' | '~*' | '!~*' | 'is' | 'is not' | 'fuzzy';

/** Value of one condition: scalar, list (`= ANY` / `!= ALL`), null (IS [NOT] NULL),
 * or dayjs (sent as ISO-8601). Always bound server-side. */
export type DomainFilterValue =
  string | number | boolean | null | Dayjs | (string | number | boolean)[];

/** One condition node. `value` may be a list (= ANY / != ALL), null (IS NULL / IS NOT NULL
 * under '='/'is'), or '@current' for the current user on user columns. Values are ALWAYS
 * bound server-side — never interpolate values into filter strings.
 *
 * The `= ANY` / `!= ALL` reading holds for columns and FK paths; on a `'<relation>.id'`
 * leaf the same shapes ask about the link SET (`!=` excludes owners linked to any of the
 * ids) — see {@link DomainRelationLink}. */
export interface DomainCondition<TColumn extends string = string> {
  /** The column tested, by name or dotted FK path. */
  property: DomainColumnRef<TColumn>;
  /** The comparison. */
  operator: DomainConditionOperator;
  /** The operand, where the operator takes one. */
  value?: DomainFilterValue;
}

/** Inside a tree, a bare string is a connector; an omitted connector means 'and'
 * (sticky — the last seen connector repeats). */
export type DomainConditionNode<TColumn extends string = string> =
  DomainCondition<TColumn> | 'and' | 'or' | DomainConditionTree<TColumn>;
/** The canonical condition tree: conditions and nested trees joined by connector strings. */
export type DomainConditionTree<TColumn extends string = string> =
  DomainConditionNode<TColumn>[];

/** Every filter slot accepts the smart-filter string (human shorthand), a single condition,
 * or the canonical condition tree — uniformly across query, queryDf, aggregate, facets,
 * count, exists, and deleteWhere (server: DomainFilterCompiler.toTree). */
export type DomainFilter<TColumn extends string = string> =
  string | DomainCondition<TColumn> | DomainConditionTree<TColumn>;

/** Runtime list of the system columns (type-level counterpart: {@link DomainSystemColumn}). */
export const DOMAIN_SYSTEM_COLUMNS = ['id', 'version', 'created_on', 'updated_on', 'author_id'] as const;

/** The per-row service columns a `withAccess` query adds (see {@link DomainQuerySpec.withAccess}):
 * booleans saying whether the CALLER may edit / delete / share that row, plus `~can_<name>` per
 * permission the table declares (`permissions` in schema.json). `~can_share` is null off row
 * mode — only row-mode tables carry per-row Share. Never exported. */
export const DOMAIN_ACCESS_COLUMNS = ['~can_edit', '~can_delete', '~can_share'] as const;

/** The per-row keys of a `withAccess` read (see {@link DOMAIN_ACCESS_COLUMNS}). */
export type DomainRowAccess = {'~can_edit': boolean; '~can_delete': boolean; '~can_share': boolean | null};

/** The soft-delete service column a `deleted: 'include' | 'only'` query projects (see
 * {@link DomainQuerySpec.deleted}): whether the SERVER holds the row as deleted. */
export const DOMAIN_DELETED_COLUMN = '~is_deleted';

/** Every service column a read can add — the {@link DOMAIN_ACCESS_COLUMNS} plus
 * {@link DOMAIN_DELETED_COLUMN}.
 * They all start with `'~'`: never exported, and hidden by `Grid.attachEditor`. */
export const DOMAIN_SERVICE_COLUMNS = [...DOMAIN_ACCESS_COLUMNS, DOMAIN_DELETED_COLUMN] as const;

/** Prefix of the caption service columns (see {@link DomainQuerySpec.captions}) — spelled once. */
export const DOMAIN_CAPTION_PREFIX = '~caption_';

/** The service column `captions: ['<column>']` projects for a ref column: the target row's
 * display name, or null where the caller cannot see the target. */
export function domainCaptionColumn(column: string): string { return `${DOMAIN_CAPTION_PREFIX}${column}`; }

/** Splits the `'<schema>.<table>'` address every domain client and UI component
 * takes, throwing on a malformed one — the single spelling of that contract. */
export function splitDomainTable(name: string): [string, string] {
  const dot = name == null ? -1 : name.indexOf('.');
  if (dot < 1 || dot === name.length - 1)
    throw new Error(`Domain table name must be '<schema>.<table>', got '${name}'`);
  return [name.substring(0, dot), name.substring(dot + 1)];
}

/** One per-column failure of a validation-failed write. */
export interface DomainColumnError { column: string; code: string; message: string; }
/** One manifest issue, addressed by the manifest path it belongs to (`'tables.orders.businessKey'`,
 * `'storage.connection'`, `'schema'`) — what {@link DomainManifestValidationError.errors}, a dry-run
 * and a live validation answer. */
export interface DomainManifestIssue {
  /** The manifest path the issue is about; absent for a schema-wide one. */
  path?: string;
  /** The server's discriminant (`'external-column-type'`, `'schema-name-taken'`, …). */
  code: string;
  /** The issue as worded for the author. */
  message: string;
}
/** Per-row failure list of a validation-failed write ({@link DomainValidationError.rows}). */
export interface DomainRowErrors { index: number; id?: string | null; errors: DomainColumnError[]; }

/** Result of one inserted row. */
export interface DomainInsertResult {
  id: string;
  created: boolean;
  version?: number;
  /** Set when the row was NOT created: 'duplicate' (business-key match) or 'idempotent-replay'. */
  status?: 'duplicate' | 'idempotent-replay';
  existingId?: string;
}
/** Result of one updated row: `version` is the new (incremented) row version. */
export interface DomainUpdateResult { id: string; version: number; }
/** Result of one deleted row. */
export interface DomainDeleteResult { id: string; deleted: true; }
/** Result of one restored row: `version` is the new (incremented) row version. */
export interface DomainRestoreResult { id: string; restored: true; version: number; }
/** Result of one transaction op (see {@link DomainOpResultFor} for the typed-tuple form). */
export type DomainOpResult = DomainInsertResult | DomainUpdateResult | DomainDeleteResult |
  DomainRestoreResult;

/** Wire shape of _audit rows (repository.dart rowAudit/tableAudit/schemaAudit selects). */
export interface DomainAuditEntry {
  /** Audit sequence number (domain_audit.id bigint — a wire number, not a uuid). */
  id: number;
  /** PG transaction id grouping multi-op writes (bigint). */
  tx_id: number | null;
  /** The row written; present on table- and schema-wide audits. */
  row_id?: string | null;
  /** The table written; present on schema audits (`'ddl'` rows carry null). */
  table_id?: string | null;
  /** What the write did. */
  op: 'insert' | 'update' | 'delete' | 'undelete' | 'promote' | 'ddl';
  /** The user who wrote; null for the system. */
  actor_id: string | null;
  /** The session the write came through. */
  session_id: string | null;
  /** Where the write came from (`'api'`, a package, a job). */
  source: string;
  /** When the write happened (ISO 8601). */
  ts: string;
  /** Column-value map before the write; null on 'insert' ops (typed `any` so existing
   * unguarded consumers compile — guard before dereferencing). */
  before: any;
  /** Column-value map after the write; null on 'delete' ops (typed `any` — see `before`). */
  after: any;
}

/** Report of domain table `deleteWhere`: rows actually soft-deleted in this call,
 * and whether more matching deletable rows remain (loop while `hasMore`). */
export interface DomainDeleteReport { deleted: number; hasMore: boolean; }

/** Report of domain table `updateWhere`: rows actually updated in this call, and whether
 * more matching editable rows remain (loop while `hasMore`). */
export interface DomainUpdateReport { updated: number; hasMore: boolean; }

/** Grantable permission on a domain registry entity (table, schema, or column schema).
 * 'Extend' is grantable on SCHEMA entities only — it lets the holder add their own
 * tables and columns to a package-managed schema the plugin opted in to. Any other
 * string is a CUSTOM permission the table declares (`permissions` in schema.json),
 * grantable on that table only and read back as `can.<name>` / `~can_<name>`. */
export type DomainPermission = 'View' | 'Edit' | 'Delete' | 'Share' | 'Extend' | (string & {});

/** One direct permission row on a domain registry entity (see `DomainTableClient.grants`). */
export interface DomainGrant {
  /** The group granted; `personal` for a user's own group. */
  group: {id: string; friendlyName: string; personal: boolean};
  /** The permission granted. */
  permission: DomainPermission;
}

/** What the caller may do with one field: `readonly` fields render as text, `editable`
 * ones take input, `immutable` ones take input on a NEW row and are read-only afterwards
 * (the key columns of an external table — settable on insert, never patched). A column the
 * caller may not SEE is absent from {@link DomainAccess.fields} altogether. */
export type FieldAccess = 'editable' | 'readonly' | 'immutable';

/** How concurrent edits are guarded on a table: `'version'` — the platform's row version
 * (`expectedVersion`); `'expected'` — the old values of the changed columns (`expected`, an
 * external table); `'none'` — the table takes no guard (a read-only table, or a binding whose
 * connector cannot compare). */
export type DomainConcurrency = 'version' | 'expected' | 'none';

/** Which part of the filter grammar a table answers: `'full'` — every operator; `'basic'` — no
 * `under`, no regex (`matches`/`!matches`), no `!like`, no datetime `!=` and no bool null test
 * (an external table). */
export type DomainFilterProfile = 'full' | 'basic';

/** What a table's {@link DomainTableClient.batch} can do beyond a plain insert. */
export interface DomainBatchSupport {
  /** `mode: 'upsert'` merges by the business key (a table declaring one, and a warehouse whose
   * connector upserts). */
  upsert: boolean;
  /** `allOrNothing: false` lands the good rows and reports the bad ones. */
  partial: boolean;
  /** `validateOnly: true` answers the dry run ({@link DomainTableClient.batch} with `validateOnly`). */
  validate: boolean;
  /** Business-key duplicates are skipped and reported (`errorOnDuplicate: false`) rather than
   * always refused. */
  skipDuplicates: boolean;
}

/** What a domain table can do AT ALL, independent of the caller — the server's one answer
 * ({@link DomainAccess.support}) for every optional behaviour a client would otherwise guess
 * from the table's shape. Orthogonal to {@link DomainAccess.can}: `can` is this caller's
 * permission, `support` is the table's storage and declaration. A client installs an optional
 * affordance only when its flag is true, and refuses by name when it is false — never probes. */
export interface DomainSupport {
  /** The system columns this table physically has, in projection order. A registration that
   * declares a subset (a platform table exposed through the `Core` schema) lists only those. */
  systemColumns: DomainSystemColumn[];
  /** The engine accepts writes ({@link DomainTableClient.insert} / {@link DomainTableClient.update} /
   * {@link DomainTableClient.batch} / {@link DomainTableClient.updateWhere}). False for a read-only,
   * entity-backed or system-subset registration — which answers `403 read-only` for every write. */
  writes: boolean;
  /** Soft delete: `deleted: 'include' | 'only'` reads and the `~is_deleted` service column. */
  deleted: boolean;
  /** {@link DomainTableClient.restore} brings a soft-deleted row back ({@link writes} AND {@link deleted}). */
  restore: boolean;
  /** Writes leave an in-transaction audit trail ({@link DomainTableClient.audit} /
   * {@link DomainTableClient.auditLog}); row watch needs it too. */
  audit: boolean;
  /** The table declares a hierarchy, so {@link DomainTableClient.pathTo} and the `under` filter
   * term answer. On a flat table `pathTo` is a {@link DomainUnsupportedError}, while `under` on a
   * column that is not a tree reads as an unknown column ({@link DomainFilterError}) — no oracle. */
  ancestors: boolean;
  /** The table carries `updated_on`, so a live client can poll for the newest write. */
  probe: boolean;
  /** The change token ({@link DomainTableClient.version}) moves: every write goes through the
   * engine. False for a registration the platform itself writes, whose token never advances —
   * poll the aggregate instead. Reading the token also needs {@link DomainAccess.can}.view. */
  version: boolean;
  /** {@link DomainTableClient.watch} subscribes to change notifications. */
  watch: boolean;
  /** {@link DomainTableClient.updateWhere} writes every row a filter matches in one call. */
  updateWhere: boolean;
  /** A query's `captions` project `~caption_<column>` with the rows; false where a client
   * resolves ref captions per row instead ({@link DomainRegistryClient.resolveNames}). */
  captions: boolean;
  /** {@link DomainsClient.transaction} — one atomic batch of inserts, updates and deletes — lands on
   * this table. */
  transaction: boolean;
  /** The guard an update carries against a concurrent edit ({@link DomainConcurrency}). */
  concurrency: DomainConcurrency;
  /** The filter grammar the table answers ({@link DomainFilterProfile}). */
  filters: DomainFilterProfile;
  /** What {@link DomainTableClient.batch} can do beyond a plain insert ({@link DomainBatchSupport}). */
  batch: DomainBatchSupport;
}

/** What {@link DomainsDataSource.draft} drafts over: an external database the caller may
 * introspect and query (`connection` is its nqName or id), one remote schema (and catalog where
 * the warehouse has them), and optionally the remote tables to draft — absent, every bindable one. */
export interface DomainDraftRequest {
  /** The DataConnection to draft over: its nqName or id. */
  connection: string;
  /** The remote schema whose tables are drafted. */
  schema: string;
  /** The remote catalog (database) holding the schema, where the warehouse has catalogs. */
  catalog?: string;
  /** The remote tables to draft; absent, every bindable table of the schema. */
  tables?: string[];
}

/** Where a manifest's rows live — its root `storage` key; absent, the platform's own store. */
export interface DomainManifestStorage {
  /** `'domain'` for the platform's own store, `'external'` for a bound warehouse. */
  kind: 'domain' | 'external';
  /** The bound DataConnection's nqName (external). */
  connection?: string;
  /** The remote schema (external). */
  schema?: string;
  /** The remote catalog (external, where the warehouse has catalogs). */
  catalog?: string;
  /** Whether rows may be inserted, updated and deleted through the connection (external). */
  writable?: boolean;
}

/** One column of a manifest table: its platform type and the keys the manifest vocabulary
 * allows (`required`, `ref`, `friendlyName`, `semType`, …). */
export interface DomainManifestColumn {
  /** The platform type (`string`, `int`, `ref`, …). */
  type?: string;
  /** The target table of a `ref` column. */
  ref?: string;
  /** Whether a row must carry a value. */
  required?: boolean;
  /** The caption shown for the column. */
  friendlyName?: string;
  /** The remote column name where it differs from the key (external). */
  column?: string;
  /** Any other key of the column vocabulary, kept as written. */
  [key: string]: unknown;
}

/** One table of a manifest, keyed by its logical name: its columns, business key and the
 * table-level keys of the manifest vocabulary. */
export interface DomainManifestTable {
  /** The columns by logical name. */
  columns: {[name: string]: DomainManifestColumn};
  /** The columns that identify a row — a bound table's primary key. */
  businessKey?: string[];
  /** The caption shown for the table. */
  friendlyName?: string;
  /** The remote table name where it differs from the logical one (external). */
  table?: string;
  /** Opted out of writes under a writable binding (external). */
  readOnly?: boolean;
  /** Any other key of the table vocabulary, kept as written. */
  [key: string]: unknown;
}

/** A domain manifest as JSON — the `domain.json` a package ships, or what {@link DomainsDataSource.draft}
 * answers and {@link DomainsDataSource.createSchema} takes; keyed by the manifest vocabulary,
 * with the root keys the platform reads named here. */
export interface DomainManifest {
  /** The schema name; the create takes it from its own argument. */
  name?: string;
  /** The manifest version, server-managed for user schemas. */
  version?: string;
  /** Where the rows live; absent, the platform's own store. */
  storage?: DomainManifestStorage;
  /** The tables by logical name. */
  tables: {[logical: string]: DomainManifestTable};
  /** Any other root key of the manifest vocabulary (`extensible`, `propertySchemas`, `migrations`, …), kept as written. */
  [key: string]: unknown;
}

/** What {@link DomainsDataSource.createSchema} takes beside the name. */
export interface DomainCreateSchemaOptions {
  /** The caption shown for the schema. */
  friendlyName?: string;
  /** A description kept on the registry entry. */
  description?: string;
  /** The manifest to register; absent, an empty platform-stored schema. */
  manifest?: DomainManifest;
  /** Every check, nothing registered. */
  dryRun?: boolean;
}

/** One remote table of the schema — the inventory lists EVERY table, requested or not, so an
 * author sees what else could be drafted. A `bindable` one is in the manifest under its
 * `logical` name (when requested, or when no `tables` were named) with the primary key as
 * `key`; the others say why not (`external-table-view`, `external-key-missing`,
 * `external-column-type`; `external-table-missing` for a requested name the schema lacks). */
export interface DomainDraftTable {
  /** The table name as the warehouse reports it. */
  remote: string;
  /** Its logical name in the manifest; absent while it is not drafted. */
  logical?: string;
  /** Whether a binding can carry it — a primary key of supported types, not a view. */
  bindable: boolean;
  /** The primary key columns, the manifest's business key. */
  key?: string[];
  /** Why it is not bindable. */
  code?: string;
  /** The reason as worded for the author. */
  message?: string;
}

/** A remote column of a drafted table with its warehouse type and, where one serves it, the platform
 * `type` the draft maps it to; `code` and `message` only on one the draft left out (a type the binding
 * cannot carry — `bytea`, `bigint` …, a malformed name, a `@param` name another column takes). */
export interface DomainDraftColumn {
  /** The remote table the column belongs to. */
  table: string;
  /** The column name as the warehouse reports it. */
  remote: string;
  /** The warehouse type (`varchar`, `int4`, …). */
  dbType: string;
  /** The platform type the draft maps it to; absent where none serves it. */
  type?: string;
  /** Why the draft left it out. */
  code?: string;
  /** The reason as worded for the author. */
  message?: string;
}

/** One foreign key touching a drafted table: `status: 'ref'` became a ref column in the manifest,
 * `'plain'` stayed a scalar and `code` says why (`external-ref-out` — the target is not in the
 * draft, `external-ref-composite`, `external-ref-key`, `external-ref-type`, `external-ref-ambiguous`). */
export interface DomainDraftRelation {
  /** The remote table holding the foreign key. */
  table: string;
  /** The foreign-key column. */
  column: string;
  /** The remote table the key points at. */
  targetTable: string;
  /** The column it points at. */
  targetColumn: string;
  /** `'ref'` where the column became a ref, `'plain'` where it stayed a scalar. */
  status: 'ref' | 'plain';
  /** Why it stayed a scalar. */
  code?: string;
  /** The reason as worded for the author. */
  message?: string;
}

/** A manifest drafted over an external database ({@link DomainsDataSource.draft}): the manifest as
 * the author's starting point, the inventory of what the warehouse has and what the draft did with
 * it, and warehouse-level diagnostics: `external-relations-unavailable` (the foreign keys could not
 * be read), `external-schema-empty` (the schema is listed but holds no table), `external-schema-missing`
 * (the connection does not list it — on PostgreSQL an empty schema reads as missing too),
 * `external-schema-unlisted` (no table, and the schema list could not be read). Nothing is
 * registered until {@link DomainsDataSource.createSchema}. */
export interface DomainDraft {
  /** The drafted manifest, the author's starting point. */
  manifest: DomainManifest;
  /** Every table, column and foreign key the warehouse has, with what the draft did with it. */
  inventory: {tables: DomainDraftTable[]; columns: DomainDraftColumn[]; relations: DomainDraftRelation[]};
  /** Warehouse-level findings, by code. */
  diagnostics: DomainManifestIssue[];
}

/** A binding's recorded verdict (ExternalBinding.STATUS_*). */
export type DomainBindingStatus = 'ok' | 'unvalidated' | 'drifted' | 'connection-missing';

/** The verdict of a schema dry run ({@link DomainsDataSource.createSchema} with `dryRun`): `'ok'`
 * with no issues — a refused manifest rejects with the typed error instead. */
export interface DomainSchemaDryRun {
  /** Always `'ok'` — a refusal rejects instead. */
  status: 'ok';
  /** Empty: the checks that would have failed reject the dry run. */
  issues: DomainManifestIssue[];
}

/** What creating a schema answers: its registry identity, physical schema (`usr_<name>` for a
 * platform-stored one, `ext_<name>` for an external binding) and, for a binding, the live
 * validation it passed. */
export interface DomainSchemaCreated {
  /** The registry id of the schema. */
  id: string;
  /** The schema name, the registry key. */
  name: string;
  /** The physical PostgreSQL schema: `usr_<name>` or `ext_<name>`. */
  pgSchema: string;
  /** The manifest version registered: `'1'` for a new schema. */
  version: string;
  /** The schema's incarnation ({@link DomainRegisteredManifest.incarnation}): new on every create of
   * the name. */
  incarnation: string;
  /** The live validation an external binding passed on create. */
  binding?: {status: DomainBindingStatus; validatedOn: string};
}

/** The recorded outcome of validating an external schema's binding against its warehouse
 * ({@link DomainSchemaClient.validate}). */
export interface DomainSchemaValidation {
  /** The schema validated. */
  schema: string;
  /** The verdict recorded on the binding. */
  status: DomainBindingStatus;
  /** When the validation ran (ISO 8601). */
  validatedOn: string;
  /** What drifted, by manifest path; empty when `status` is `'ok'`. */
  issues: DomainManifestIssue[];
}

/** A registered schema's manifest as {@link DomainSchemaClient.manifest} answers it: the manifest
 * vocabulary reconstructed from the registry, plus the two tokens an apply sends back as
 * `ifVersion` / `ifIncarnation`. */
export interface DomainRegisteredManifest extends DomainManifest {
  /** The schema name. */
  name: string;
  /** The apply counter: what {@link DomainApplyBody.ifVersion} names. */
  version: string;
  /** The schema row's creation instant (ISO 8601): what {@link DomainApplyBody.ifIncarnation} names.
   * A name deleted and re-created starts at version `'1'` again with a new incarnation. */
  incarnation: string;
}

/** A group as the access surfaces name it. */
export interface DomainGroupRef {
  /** The group id — what a grant is addressed to. */
  id: string;
  /** The group's caption; a user's name for a personal group. */
  friendlyName: string;
  /** Whether it is a user's own group. */
  personal: boolean;
}

/** One group's COMPLETE direct permission set on a registry entity, as the access snapshot
 * lists it — Share and custom permissions included, nothing inherited. */
export interface DomainGrantSet {
  /** The group granted. */
  group: DomainGroupRef;
  /** Every permission the group holds directly on the entity. */
  permissions: DomainPermission[];
}

/** The direct grants on one registry entity in the access snapshot. Answered only where the
 * caller holds Share on it — elsewhere `canShare` is false and `grants` is `null`: unknown, not
 * empty. */
export interface DomainEntityAccess {
  /** The registry entity id. */
  entityId: string;
  /** Whether the caller may change the entity's grants (Share on it). */
  canShare: boolean;
  /** The direct grants, complete per group; `null` when the caller may not read them. */
  grants: DomainGrantSet[] | null;
}

/** One table of the access snapshot: its own grants and its core schema's (the
 * everyone-visible property schema its unrestricted columns belong to — informational: an apply
 * edits table grants and column restrictions only). */
export interface DomainTableAccess extends DomainEntityAccess {
  /** The physical (remote, for a binding) table name. */
  remote: string;
  /** The table's core property schema — the Share target of its column restrictions. */
  coreSchema: DomainEntityAccess;
}

/** One relational column of the access snapshot: whether it is split out of the core schema and,
 * when it is, who reads and edits it. The lists are `null` (unknown) where the caller may not
 * Share the table's core schema. */
export interface DomainColumnAccess {
  /** `'restricted'` when the column has its own per-column schema. */
  state: 'unrestricted' | 'restricted';
  /** Whether the caller may restrict or unrestrict the column (Share on the core schema). */
  canShare: boolean;
  /** The per-column schema id of a restricted column the caller may Share. */
  schemaId?: string;
  /** The ACL revision of a restricted column the caller may Share: it moves with every grant,
   * revoke and re-restriction, so an unrestriction sent with it
   * ({@link DomainColumnUnrestriction.revision}) cannot undo a change made after the snapshot. */
  revision?: string;
  /** Groups holding View on a restricted column. */
  view?: DomainGroupRef[] | null;
  /** Groups holding Edit on a restricted column. */
  edit?: DomainGroupRef[] | null;
  /** Any other permission held on a restricted column, per group. */
  other?: {group: DomainGroupRef; permission: DomainPermission}[] | null;
}

/** The lossless access snapshot of a schema for its editor ({@link DomainSchemaClient.access}),
 * read in one transaction with the tokens it stands for: per table the direct grants, per column
 * the restriction state. `null` anywhere means the caller may not read that part (no Share on the
 * target), never that it is empty. */
export interface DomainSchemaAccess {
  /** Per table, by logical name. */
  tables: {[logical: string]: DomainTableAccess};
  /** Per relational column, by `'<table>.<column>'`. */
  columns: {[address: string]: DomainColumnAccess};
  /** The apply counter the snapshot was read at ({@link DomainApplyBody.ifVersion}). */
  version: string;
  /** The schema's incarnation ({@link DomainApplyBody.ifIncarnation}). */
  incarnation: string;
}

/** One permission of one group on one table — the unit of {@link DomainAccessDelta.grant} and
 * {@link DomainAccessDelta.revoke}. */
export interface DomainGrantTriple {
  /** The table's logical name. */
  table: string;
  /** The group id. */
  group: string;
  /** The permission granted or revoked. */
  permission: DomainPermission;
}

/** One permission of one group, inside a column restriction. */
export interface DomainGroupPermission {
  /** The group id. */
  group: string;
  /** `'View'` or `'Edit'` — what a column restriction grants. */
  permission: DomainPermission;
}

/** A column addressed by table and column name. */
export interface DomainColumnTarget {
  /** The table's logical name. */
  table: string;
  /** The column's logical name. */
  column: string;
}

/** Restricts a column (splits it out of the everyone-visible core schema, idempotently) and
 * applies grant/revoke DELTAS on its per-column schema — never a full list, so a grant somebody
 * else made meanwhile survives an edit that never named it. A key column of an external table
 * cannot be restricted (`'key-column'`). */
export interface DomainColumnRestriction extends DomainColumnTarget {
  /** The restriction state the edit was made from. A standalone unrestriction moves no schema
   * version, so this is what catches it: sent, the apply is refused under the lock
   * (`'access-conflict'`, naming the column and both states) when the column's state differs;
   * omitted, the op applies whatever the state. */
  from?: 'restricted' | 'unrestricted';
  /** The {@link DomainColumnAccess.revision} the edit was made from; sent, the apply is refused
   * (`'access-conflict'`, with `expectedRevision` and `currentRevision`) when the column's grants
   * changed since. */
  revision?: string;
  /** Permissions to grant on the per-column schema. */
  grant?: DomainGroupPermission[];
  /** Permissions to revoke from the per-column schema. */
  revoke?: DomainGroupPermission[];
}

/** Makes a column visible to everyone again: its per-column schema and the grants on it go. */
export interface DomainColumnUnrestriction extends DomainColumnTarget {
  /** `'restricted'`, the state the edit was made from; refused (`'access-conflict'`) when the
   * column is not restricted any more. Omitted, the op is idempotent. */
  from?: 'restricted';
  /** The {@link DomainColumnAccess.revision} the edit was made from: refused (`'access-conflict'`,
   * with `expectedRevision` and `currentRevision`) when anyone granted, revoked, or re-restricted
   * the column since — an unrestriction is not a delta, so a stale one would otherwise restore
   * access revoked meanwhile. */
  revision?: string;
}

/** The access changes an apply carries ({@link DomainApplyBody.access}): exact permission-triple
 * deltas, executed in the registry transaction under the deploy lock AFTER the manifest change,
 * so a column added and restricted in one apply is never visible to today's readers and
 * `writable` turned on never activates a grant the same apply revokes. Every op needs Share on
 * its target (a column's: the table's core schema), or the whole apply is refused
 * (`'access-forbidden'`, naming the targets); a target the same apply drops is refused too
 * (`'access-target-dropped'` — the purge owns it). A triple both granted and revoked, or a column
 * named by two restriction ops, is `'invalid-access'`. */
export interface DomainAccessDelta {
  /** Table grants to add. */
  grant?: DomainGrantTriple[];
  /** Table grants to remove. */
  revoke?: DomainGrantTriple[];
  /** Columns to restrict, with the deltas on their readers and editors. */
  restrict?: DomainColumnRestriction[];
  /** Columns to make visible to everyone again (their per-column schema and its grants go). */
  unrestrict?: DomainColumnUnrestriction[];
}

/** What {@link DomainSchemaClient.apply} takes: a partial manifest — named tables replace their
 * registered definition WHOLESALE, everything omitted stays as registered — with the two tokens
 * the edit was made against, the schema-level metadata, and the access deltas. */
export interface DomainApplyBody {
  /** The apply counter the edit was made against ({@link DomainRegisteredManifest.version}):
   * REQUIRED on a user-managed schema (`'version-required'` without it), checked under the deploy
   * lock — of two concurrent applies exactly one commits, the other rejects with a
   * {@link DomainVersionConflictError}. On a package-managed schema (the extension path) it is
   * optional and tracks `ext_version`. */
  ifVersion?: string;
  /** The incarnation the edit was made against ({@link DomainRegisteredManifest.incarnation}): a
   * token from another life of the name (deleted and re-created since) rejects with a
   * {@link DomainVersionConflictError} carrying both incarnations. Optional on the wire; send it. */
  ifIncarnation?: string;
  /** Tables to create or replace, whole descriptors by logical name. */
  tables?: {[logical: string]: DomainManifestTable};
  /** Tables to unregister (a destructive plan: needs `confirmDestructive`, and Delete on each). */
  dropTables?: string[];
  /** Package-managed schemas only: your own columns on plugin tables that opted in, as the full
   * state of YOUR columns of each table. */
  extend?: {[table: string]: {columns: {[column: string]: DomainManifestColumn}}};
  /** Property schemas to merge by name. */
  propertySchemas?: {[name: string]: unknown};
  /** The schema's caption (user-managed schemas). */
  friendlyName?: string;
  /** The schema's description (user-managed schemas). */
  description?: string;
  /** `writable` is the one storage key an apply may change, and only on an external binding;
   * any other key must repeat the registered value (`'storage-immutable'`; the kind cannot flip,
   * `'storage-conversion'`). */
  storage?: {writable?: boolean};
  /** Access changes executed with the manifest change, atomically. */
  access?: DomainAccessDelta;
  /** Required when the plan is destructive (`'destructive-confirmation-required'` otherwise, the
   * plan riding in the error's body). */
  confirmDestructive?: boolean;
}

/** The effect one access op had, or would have: `'none'` where the triple was already there, or
 * already gone. */
export interface DomainGrantEffect {
  /** The group, resolved. */
  group: DomainGroupRef;
  /** The permission the op named. */
  permission: DomainPermission;
  /** What the op does to the permission row. */
  effect: 'grant' | 'revoke' | 'none';
}

/** The effect of a table grant or revoke. */
export interface DomainTableGrantEffect extends DomainGrantEffect {
  /** The table's logical name. */
  table: string;
}

/** The effect of a column restriction: the split itself, then each grant and revoke on the
 * per-column schema. */
export interface DomainRestrictEffect extends DomainColumnTarget {
  /** `'restrict'` when the column is split out by this apply; `'none'` when it already was. */
  effect: 'restrict' | 'none';
  /** The effect of each grant named. */
  grants: DomainGrantEffect[];
  /** The effect of each revoke named. */
  revokes: DomainGrantEffect[];
}

/** The effect of making a column visible to everyone again. */
export interface DomainUnrestrictEffect extends DomainColumnTarget {
  /** `'unrestrict'` when the per-column schema goes; `'none'` when there was none. */
  effect: 'unrestrict' | 'none';
}

/** The access section of a plan ({@link DomainApplyPlan.access}): every op with its effect — in a
 * dry run as planned against the rows as they are now plus what the apply bootstraps, in an
 * applied answer as resolved under the lock. */
export interface DomainAccessEffects {
  /** The table grants. */
  grant: DomainTableGrantEffect[];
  /** The table revokes. */
  revoke: DomainTableGrantEffect[];
  /** The column restrictions. */
  restrict: DomainRestrictEffect[];
  /** The columns made visible again. */
  unrestrict: DomainUnrestrictEffect[];
}

/** What the purge of a dropped table takes with it, as counts — names would leak what the
 * reviewer may not see. */
export interface DomainLostTable {
  /** Direct grants on the table. */
  grants: number;
  /** Grants on the table's core property schema. */
  coreSchemaGrants: number;
  /** Per restricted column, the grants on its per-column schema. */
  restrictions: {[column: string]: number};
  /** Rows promoted to entities (individually shareable). */
  promotedRows: number;
  /** Grants on those promoted rows. */
  rowGrants: number;
  /** Saved filters of the table — DELETED with it. */
  savedFilters: number;
  /** Saved filters of other tables that travel into this one through a ref — kept, but no
   * longer resolving. */
  affectedFilters: number;
}

/** What dropping a column takes with it. */
export interface DomainLostColumn {
  /** Whether the column had its own per-column schema. */
  restricted: boolean;
  /** Grants on that per-column schema. */
  grants: number;
  /** Saved filters naming the column (by column, by FK path, inside an expression) — kept, but
   * no longer resolving. */
  affectedFilters: number;
}

/** What a user-schema apply loses ({@link DomainApplyPlan.lost}): per dropped table, per dropped
 * column, and per re-pointed or demoted ref (`'<table>.<column>'`) the saved filters travelling
 * through it, however many refs the path crosses (up to four). Each entry counts a filter once, but
 * the entries overlap: a filter travelling through a demoted ref into a dropped table counts under
 * both `refs` and `tables`. */
export interface DomainApplyLost {
  /** Per dropped table, by logical name. */
  tables: {[table: string]: DomainLostTable};
  /** Per dropped column, by `'<table>.<column>'`. */
  columns: {[address: string]: DomainLostColumn};
  /** Per ref column whose type changed, by `'<table>.<column>'`; `demoted` when it becomes a plain
   * value — its own values stay, only the paths through it break. */
  refs: {[address: string]: {affectedFilters: number; demoted?: true}};
}

/** The change plan of {@link DomainSchemaClient.apply}: what a dry run answers, and what a commit
 * answers with `applied: true` (and the access effects as resolved under the lock). On an
 * external binding the plan is `registrationOnly` — registry rows move, the warehouse is never
 * touched, no live counts. */
export interface DomainApplyPlan {
  /** The manifest version the apply registers. */
  version: string;
  /** Whether the plan drops data, constraints or registrations — `confirmDestructive` needed. */
  destructive: boolean;
  /** What the apply creates: tables by name, and per table its new columns, uniques and constraints. */
  creates: {tables: string[]; columns: {[table: string]: string[]}; uniques: {[table: string]: string[]};
    constraints: {[table: string]: string[]}};
  /** What the apply drops; `liveRows` / `nonNullValues` on a platform-stored schema only. */
  drops: {tables: {name: string; liveRows?: number}[]; columns: {table: string; column: string; nonNullValues?: number}[]};
  /** Columns whose type changes. */
  typeChanges: {table: string; column: string; to: string}[];
  /** Unique constraints dropped. */
  uniqueDrops: {table: string; column: string}[];
  /** Columns turned required or optional, with the rows a NOT NULL would refuse. */
  requiredToggles: {table: string; column: string; required: boolean; violatingRows: number}[];
  /** Tables whose business key changes. */
  businessKeyChanges: {table: string; businessKey: string[]}[];
  /** Per table, the columns whose autonumbering changes. */
  autoNumberChanges: {[table: string]: string[]};
  /** Table-level toggles (audit, name column, relations, permissions), per kind. */
  alters: {[kind: string]: unknown};
  /** Declared constraints dropped. */
  constraintDrops: unknown;
  /** Data pre-checks a destructive change would violate. */
  violations: unknown[];
  /** Manifest issues that refuse the plan (empty on an answered plan — a refusal rejects). */
  refusals: DomainManifestIssue[];
  /** Restricted columns of an external target that are now key columns — visible to everyone
   * again, since the row id encodes them. */
  keyColumnsUnrestricted: DomainColumnTarget[];
  /** The registered-snapshot transition with its migration scaffold; null when none applies. */
  migration: {from: string; to: string; changes: unknown[]; up: string; down: string} | null;
  /** External binding: registry rows only, no warehouse DDL and no row counts. */
  registrationOnly?: true;
  /** User-managed schema: what the drops purge or break. */
  lost?: DomainApplyLost;
  /** User-managed schema: the schema-level metadata that changes, `{from, to}` per key
   * (`friendlyName`, `description`), so a rename-only plan is not empty. */
  metadata?: {[key: string]: {from: string | null; to: string | null}};
  /** The access effects, present when the body carried `access`. */
  access?: DomainAccessEffects;
  /** Set on the answer of a commit; absent on a dry run. `false` with {@link noop}: the apply
   * changed nothing, so nothing was written. */
  applied?: boolean;
  /** Set on the answer of a commit that found nothing to change — the schema already holds what
   * the body declares: nothing is written and `version` stays the one the apply was made against. */
  noop?: true;
  /** Package-managed schema: the extension version the apply moved to. */
  extVersion?: number;
}

/** The answer of a committed {@link DomainSchemaClient.apply}: the plan, `applied`, with the
 * access effects as resolved under the lock — or, when the apply changed nothing, `applied: false`
 * with `noop: true` and the version unchanged. */
export type DomainApplied = DomainApplyPlan & ({applied: true} | {applied: false; noop: true});

/** Effective access of the CURRENT user on one domain table
 * (see `DomainTableClient.access`). Composed by the SERVER
 * (`GET /domains/{schema}/{table}/access`) from the same predicates its reads
 * and writes apply: grants on the FINAL securing entity through the master delegate
 * chain, column security, and relation travel.
 *
 * `can` holds TABLE-level affordances, and they drift toward denial: every flag is
 * derived from grants on the securing entity, so grants that reach individual ROWS
 * (a master row of a master-mode table, a promoted row of a row-mode table) are not
 * counted. A false flag therefore means "no table-wide right" — the caller may still
 * succeed on a particular row, and per-row truth rides the rows of a
 * `withAccess` query ({@link DOMAIN_ACCESS_COLUMNS}). A true flag mirrors a predicate
 * the server also enforces, so gating UI on it avoids the common 403s (but never
 * assume it removes them: grants can change between the probe and the write, and
 * column-level restrictions are checked per value). Read-only registrations
 * (platform tables exposed through the `Core` schema) answer every write flag false
 * and every field `readonly` for everyone, admins included. */
export interface DomainAccess {
  /** What the caller may do, capability by capability. */
  can: {
    /** View grant on the securing entity. Row-mode tables may still expose
     * individually granted rows when false. */
    view: boolean;
    /** Edit grant on the securing entity — the server's insert predicate — on a
     * table that accepts writes. False negative on master-mode tables where the
     * caller holds the grant on an individual MASTER ROW rather than the master table. */
    insert: boolean;
    /** {@link insert} AND at least one column is editable for the caller (the
     * built-in grid-editability rule) AND the table grant actually reaches rows —
     * for non-table security modes that needs the securing table's rows to default
     * to table visibility, otherwise access is per-row. Same false-negative shape.
     * Computed before the `immutable` rewrite of {@link FieldAccess}: on an external
     * table whose only writable columns are its keys it is true although no
     * saved row has an editable field. */
    edit: boolean;
    /** Delete grant on the securing entity, under the same reaches-rows rule as
     * {@link edit}; same false-negative shape. */
    delete: boolean;
    share: boolean;
    /** A permission the table declares (`permissions` in schema.json), by name:
     * a grant of it on the securing entity, under the writability rules above. */
    [custom: string]: boolean;
  };
  /** Every column the caller may see (declared and system columns alike), keyed by
   * name; `editable` iff the caller may write it (Edit on an owning property
   * schema); `immutable` for a key column of an external table the caller may write
   * (settable on insert, never patched). A restricted column is ABSENT. An `autoNumber` column reads `readonly`
   * although the server still accepts a supplied value (imports keep their numbers)
   * — a form never sends one. */
  fields: {[column: string]: FieldAccess};
  /** Relations the caller may expand (View on both the junction and the target
   * table), in declaration order. */
  travelableRelations: string[];
  /** The FINAL securing table (`<schema>.<table>`) every grant above is evaluated on:
   * the table itself, or the end of its master delegate chain. */
  securingTable: string;
  /** How rows are secured: by the table's own grants, a master row's, or each row's. */
  securityMode: 'table' | 'master' | 'row';
  /** Whether writes leave an in-transaction audit trail. */
  audit: boolean;
  /** Whether the table declares a natural key, which upserts match on. */
  hasBusinessKey: boolean;
  /** What the TABLE can do, independent of this caller — computed once on the server from the
   * registration, so a client never guesses an affordance from the table's shape. */
  support: DomainSupport;
}

/** FK-inverted reference to a child (detail) table — drives detail-table links
 * and entity-view child tabs (see `DomainRegistryClient.tableInfo`). */
export interface DomainChildTableRef {
  /** The child table's schema. */
  schema: string;
  /** The child table. */
  table: string;
  /** The child table's FK column referencing this table. */
  fkColumn: string;
  /** Friendly label of the FK column ('issue_id' → 'Issue') — disambiguates
   * children referencing the same parent through several columns. */
  label: string;
}

/** Registry metadata of one domain table (see `grok.dapi.domains.registry`). */
export interface DomainTableInfo {
  /** Primary display-name column; null when the table declares none. */
  nameColumn: string | null;
  /** Natural-key columns; empty when the table declares none. */
  businessKey: string[];
  /** How an app spells a row in a URL: `'businessKey'` — the business-key values (a platform
   * table); `'id'` — the canonical row id as one opaque segment (an external table, whose key
   * values may carry any delimiter). */
  rowAddress: 'id' | 'businessKey';
  /** The DECLARED `friendlyName` (`friendlyName` in schema.json); undefined when
   * the table declares none — unlike {@link singularName}/{@link pluralName},
   * nothing is derived here, so the caller can tell a chosen label apart from
   * the capitalized table name. */
  friendlyName?: string;
  /** Effective singular display name (declared, else derived from the table name). */
  singularName: string;
  /** Effective plural display name (declared, else derived from the table name). */
  pluralName: string;
  /** How rows are secured: by the table's own grants, a master row's, or each row's. */
  securityMode: 'table' | 'master' | 'row';
  /** Whether writes leave an in-transaction audit trail. */
  audit: boolean;
  /** The table declares a tree (`hierarchy` in schema.json): exactly one ref column
   * targets the table itself, and {@link DomainTableClient.pathTo} / the `under`
   * filter term walk it. */
  hierarchy: boolean;
  /** The self-referencing column the tree is built on; null off a hierarchy table. */
  parentColumn: string | null;
  /** What the table holds, as declared in schema.json; null when undeclared. */
  description: string | null;
  /** The tables whose foreign keys point at this one. */
  childTables: DomainChildTableRef[];
  /** The columns {@link DomainQuerySpec.search} matches: those declared
   * `searchable: true`, else the name column, else none. */
  searchableColumns: string[];
  /** Declared table constraints (`constraints` in schema.json), in declaration order. */
  constraints: {name: string; expr: string; message?: string}[];
  /** Per ref column, the declared `filter` its candidates are narrowed by (a
   * smart-filter string; `$<column>` binds a sibling value of the row being edited). */
  refFilters: {[column: string]: string};
  /** Custom permission names the table declares (`permissions` in schema.json). */
  permissions: string[];
}

/** Per-row outcome inside a {@link DomainBatchReport}. */
export interface DomainBatchRowResult {
  /** The row's position in the batch. */
  index: number;
  /** The row's id; null where nothing landed. */
  id: string | null;
  /** `'merged'`: an upsert landed on a storage that cannot tell an insert from an update. */
  status: 'inserted' | 'updated' | 'merged' | 'duplicate' | 'error';
  /** The id of the row a `'duplicate'` matched. */
  existingId?: string;
  /** The per-column failures of an `'error'` row. */
  errors?: DomainColumnError[];
}

/** Base of all structured /domains failures; [body] is the server's error envelope. */
export class DomainError extends Error {
  constructor(message: string, public readonly status: number,
              public readonly body: {[key: string]: any}) {
    super(message);
    this.name = new.target.name;
  }
  /** Server discriminant: 'validation' | 'version-conflict' | 'access-conflict' | 'restrict' | 'filter' |
   * 'forbidden' | 'not-found' | 'unsupported' | 'manifest-validation' | 'invalid-mode' | 'id-collision' |
   * 'destructive-confirmation-required' | 'not-user-managed' | 'not-package-managed' |
   * 'schema-management-unavailable' | 'apply-failed' | 'delete-failed' | '' (transport). */
  get code(): string { return `${this.body['error'] ?? ''}`; }
  /** Index of the failing op in a transaction ops list; undefined otherwise. */
  get opIndex(): number | undefined { return this.body['opIndex']; }
}

export class DomainValidationError extends DomainError {          // code 'validation'
  get rows(): DomainRowErrors[] { return this.body['rows'] ?? []; }
  /** True when the failure is a business-key/idempotency uniqueness conflict
   * (per-row error code 'unique' — how the server reports duplicates under errorOnDuplicate). */
  get isDuplicate(): boolean {
    return this.rows.some((r) => (r.errors ?? []).some((e) => e.code === 'unique'));
  }
}
export class DomainVersionConflictError extends DomainError {     // code 'version-conflict'
  get id(): string { return this.body['id']; }
  /** The row's version as the server holds it; a schema apply's conflict carries the schema's
   * apply counter, a string. */
  get currentVersion(): number | string { return this.body['currentVersion']; }
  /** The version the write named. */
  get expectedVersion(): number | string { return this.body['expectedVersion']; }
  /** The `ifIncarnation` a schema apply sent; present whenever the apply sent one, undefined on a
   * row conflict. It differs from {@link currentIncarnation} only when the name was deleted and
   * re-created since the edit was made. */
  get expectedIncarnation(): string | undefined { return this.body['expectedIncarnation']; }
  /** The schema's current incarnation, beside {@link expectedIncarnation}. */
  get currentIncarnation(): string | undefined { return this.body['currentIncarnation']; }
  /** The old-value guard the write carried ({@link DomainTransactionOp.expected}), per column;
   * undefined when the conflict is a version one. */
  get expected(): {[column: string]: unknown} | undefined { return this.body['expected']; }
  /** The guard columns as the server read them AFTER refusing the write — a later observation,
   * they may have moved again. */
  get current(): {[column: string]: unknown} | undefined { return this.body['current']; }
}
export class DomainRestrictError extends DomainError {}           // code 'restrict'
export class DomainFilterError extends DomainError {}             // code 'filter'
export class DomainForbiddenError extends DomainError {}          // code 'forbidden'
export class DomainNotFoundError extends DomainError {}           // code 'not-found'
export class DomainManifestValidationError extends DomainError {  // code 'manifest-validation'
  /** The issues by manifest path ({@link DomainManifestIssue}). */
  get errors(): DomainManifestIssue[] { return this.body['errors'] ?? []; }
}
/** The one refusal for an operation this table's storage or declaration cannot do (422) —
 * `pathTo` on a flat table. Distinct from a permission failure
 * ({@link DomainForbiddenError}): no grant makes it succeed. {@link DomainAccess.support}
 * says the same thing BEFORE the call, so a client should never see this one. */
export class DomainUnsupportedError extends DomainError {         // code 'unsupported'
  /** The operation the table cannot do ('ancestors', 'batch', …). */
  get op(): string { return `${this.body['op'] ?? ''}`; }
}

/** Maps the interop's '#domainError' plain object to a typed class; passes anything else through. */
export function toDomainError(e: any): any {
  if (e == null || e['#domainError'] !== true)
    return e;
  const body = e.body ?? {};
  const ctor = ({
    'validation': DomainValidationError,
    'version-conflict': DomainVersionConflictError,
    'restrict': DomainRestrictError,
    'filter': DomainFilterError,
    'forbidden': DomainForbiddenError,
    'not-found': DomainNotFoundError,
    'unsupported': DomainUnsupportedError,
    'manifest-validation': DomainManifestValidationError,
  } as {[code: string]: typeof DomainError})[`${body['error']}`] ?? DomainError;
  return new ctor(`${e.message ?? body['message'] ?? 'Domain request failed'}`, e.status ?? 0, body);
}

/** Awaits [p], rethrowing interop-marked failures as typed DomainErrors. */
export async function domainCall<T>(p: Promise<T>): Promise<T> {
  try {
    return await p;
  } catch (e: any) {
    throw toDomainError(e);
  }
}

/** Parameter values of the `DomainQuery` function — the serializable representation of
 * what a domain view is showing (filter elements, ordering, cap), as
 * {@link DomainView.query} reports it. Feeding it back to the function
 * ({@link DomainQuery.run}) reproduces the subset as a DataFrame, with the creation
 * script recorded.
 *
 * This is the ONLY parameter shape of the function — {@link DomainQuery} is its
 * mutable, URL-serializable counterpart. Grammar strings, not condition trees: each
 * element is one `DomainQuery` filter expression (`'status = "open"'`) with the
 * per-element JSON escape hatch. Use {@link DomainQuerySpec} for the REST surface
 * instead.
 *
 * Aggregate mode (non-empty `aggregations`/`groupBy`) REJECTS `columns` and `offset`
 * with a 400 — they are select-mode parameters, not ignored ones. */
export interface DomainQueryParams {
  /** The queried table's schema. */
  schema: string;
  /** The queried table. */
  table: string;
  /** Projection; omit for all viewable columns (select mode only). */
  columns?: string[];
  /** Filter elements, AND-combined. OMITTED, not `[]`, when there are none —
   * these are the function's parameter values, and it takes no empty lists
   * (`DomainView.query` on an unfiltered view has no `filters` key). */
  filters?: string[];
  /** Master-FK expands (`'<fk_column>'`); `'details:'` child arrays are not supported. */
  joins?: string[];
  /** Measures (`'count'`, `'avg(amount) as avg_amount'`) — non-empty means aggregate mode. */
  aggregations?: string[];
  /** Grouping columns — non-empty means aggregate mode. */
  groupBy?: string[];
  /** Sort elements; a `'!'` prefix means descending. */
  orderBy?: string[];
  limit?: number;
  /** Row offset (select mode only). */
  offset?: number;
}

/** The part of a read that selects ROWS rather than shapes the page: what
 * {@link DomainTableClient.count}, {@link DomainTableClient.exists},
 * {@link DomainTableClient.aggregate} and {@link DomainTableClient.query} all narrow by.
 * One object, so a caller cannot pass `filter` and forget `deleted` — a total that disagrees
 * with the rows it counts is the bug this shape removes. */
export interface DomainReadScope<TColumn extends string = string> {
  /** Smart-filter string (same grammar as entity search, e.g. `barcode starts "P-1"`),
   * a single condition, or the canonical condition tree — values are bound server-side.
   * Row caps: JSON `query` 10k, d42 `queryDf` 10M. */
  filter?: DomainFilter<TColumn>;
  /** Case-insensitive substring over the table's searchable columns — those declared
   * `searchable: true`, else the name column ({@link DomainTableInfo.searchableColumns});
   * a table with neither rejects with a filter error. ANDed with `filter`. */
  search?: string;
  /** Soft-deleted rows: `'exclude'` (the default — only live rows), `'include'` (both) or
   * `'only'` (the trash). Anything but `'exclude'` also projects `~is_deleted`
   * ({@link DOMAIN_SERVICE_COLUMNS}); a deleted row is read-only until
   * {@link DomainTableClient.restore} brings it back. Needs {@link DomainSupport.deleted}. */
  deleted?: 'exclude' | 'include' | 'only';
}

/** Query options for domain table `query`/`queryDf`: a {@link DomainReadScope} plus what
 * shapes the page. */
export interface DomainQuerySpec<TColumn extends string = string, TExpandKey extends string = string>
    extends DomainReadScope<TColumn> {
  /** Comma-separated column list; `!` prefix for descending, e.g. `'name,!created_on'`. */
  sort?: string;
  /** Columns to return; omit for all viewable columns. */
  columns?: TColumn[];
  /** Expansions (same-schema, depth 1): `'<fk_column>'` returns the master row's declared
   * columns prefixed `'<fk_column>.<name>'`; `'details:<table>[.<fk_column>]'` returns capped
   * child-row arrays under the child-table name (JSON queries only — not supported by
   * `queryDf`); `'<relation>'` (a declared many-to-many) returns the link array — see
   * {@link DomainRelationLink}. */
  expand?: TExpandKey[];
  /** The row cap. */
  limit?: number;
  /** Rows skipped before the first returned. */
  offset?: number;
  /** Adds the per-row {@link DOMAIN_ACCESS_COLUMNS} (`~can_edit`, `~can_delete`, `~can_share`,
   * plus `~can_<name>` per declared custom permission): JSON rows carry them as boolean keys,
   * `queryDf` frames as trailing bool columns. Row-level truth where the table-level
   * {@link DomainAccess.can} flags are false negatives. */
  withAccess?: boolean;
  /** Ref columns to project the TARGET ROW'S DISPLAY NAME for, each as one service column
   * `'~caption_<column>'` ({@link domainCaptionColumn}) — string, nullable, in request order.
   * A list renders `project_id` as 'Apollo' without a second round trip and without the
   * caller knowing which column of the target is its name.
   *
   * Name only: never a target field (use {@link expand} for those), never nested
   * (`'a.b'` rejects), never a duplicate, and never a non-ref or invisible column — all four
   * reject with a {@link DomainFilterError} that names nothing else (no oracle). A target row
   * the caller may not View reads NULL, exactly like a target that does not exist. Captions are
   * not fields: absent from {@link DomainAccess.fields}, never editable, excluded from csv and
   * binary export. Independent of {@link columns}: a caption may be asked for a ref column that
   * is not projected, and asking for one never projects the ref column itself.
   *
   * One LEFT JOIN per caption per page — never default them on. */
  captions?: TColumn[];
}

/** One measure of {@link DomainAggregateSpec}: `fn` over `column` (`count` needs no column);
 * the output name is `as` (defaults to `fn` or `<fn>_<column>`). */
export interface DomainAggregateMeasure<TColumn extends string = string, TAlias extends string = string> {
  /** The aggregate function. */
  fn: 'count' | 'sum' | 'avg' | 'min' | 'max';
  /** The column aggregated; `count` needs none. */
  column?: DomainColumnRef<TColumn>;
  /** The output name; `fn` or `fn_column` by default. */
  as?: TAlias;
}

/** Spec for domain table `aggregate`: a {@link DomainReadScope} (the rows aggregated —
 * `search` and `deleted` narrow it exactly as they narrow {@link DomainQuerySpec}) plus the
 * grouping. */
export interface DomainAggregateSpec<TColumn extends string = string,
    TGroup extends string = string, TAlias extends string = string>
    extends DomainReadScope<TColumn> {
  /** Column names to group by; omit for a single-row grand total. */
  groupBy?: (DomainColumnRef<TColumn> & TGroup)[];
  /** The measures computed per group. */
  measures: DomainAggregateMeasure<TColumn, TAlias>[];
  /** Comma-separated output names (group columns or measure aliases); `!` prefix for descending. */
  sort?: string;
  /** The group-row cap. */
  limit?: number;
  /** Master-FK expands ('<fk_column>'); groupBy/measures may then use '<fk>.<col>' (§13.2). */
  expand?: string[];
}
/** One result row of `aggregate` — keys are the group columns and measure aliases. */
export type DomainAggregateRow<TKeys extends string = string> =
  {[K in TKeys]: number | string | boolean | null};

/** The facet kinds of {@link DomainFacetSpec}. */
export type DomainFacetKind = 'categories' | 'histogram' | 'minMax' | 'count' | 'plan';

/** One facet request of {@link DomainFacetsSpec}; `id` keys its result in the response. */
export interface DomainFacetSpec<TColumn extends string = string,
    TId extends string = string, TKind extends DomainFacetKind = DomainFacetKind> {
  /** Response key for this facet's result. */
  id: TId;
  /** The facet kind, which decides the result shape ({@link DomainFacetResultOf}). */
  kind: TKind;
  /** Column name; `'categories'` also accepts a dotted FK path (e.g. `'category_id.name'`, up to
   * 3 hops — the counts respect the referenced table's row predicate) and a declared
   * many-to-many path (`'labels.name'`, or `'labels.id'` for the id + display-name form a
   * checkbox list wants): those count DISTINCT OWNERS, so an owner with two labels counts once
   * under each, and only owners with no VISIBLE link fall in the null bucket.
   * Not used by `'count'`/`'plan'`. */
  column?: DomainColumnRef<TColumn>;
  /** `'plan'` only: columns to profile (capped distinct count plus numeric/datetime min/max). */
  columns?: TColumn[];
  /** `'categories'` only: category cap (default 100, clamped to 1..1000); the result carries
   * `hasMore` when more remain — narrow with `search` instead of raising the cap. */
  limit?: number;
  /** `'categories'` only: server-side substring narrowing (compiled as a bound ILIKE). */
  search?: string;
  /** `'histogram'` only: bucket count (default 20, clamped to 1..200). */
  bins?: number;
  /** `'histogram'` only: pinned lower bound (a number, or an ISO-8601 string for datetime columns)
   * so zooming does not re-derive the axis; omit to bucket over the data bounds. */
  min?: number | string;
  /** `'histogram'` only: pinned upper bound (see `min`). */
  max?: number | string;
}

/** Spec for domain table `facets`: one request computes every facet in a single round
 * trip. Category counts, histogram buckets, and `'count'` are computed under `filter` with the
 * conditions on that facet's own column stripped, so a filter control shows counts under all
 * OTHER filters (classic faceted search). Exception: `'minMax'`, `'plan'`, and histogram BOUNDS
 * are computed under the row predicate only — the stable-axis rule ignores `filter` so a
 * narrowing filter never re-derives the axis. At most 32 facets per request. */
export interface DomainFacetsSpec<TColumn extends string = string,
    TId extends string = string, TKind extends DomainFacetKind = DomainFacetKind> {
  /** Smart string, single condition, or condition tree; omit for unfiltered counts. */
  filter?: DomainFilter<TColumn>;
  /** The facets computed in one round, each answered under its `id`. */
  facets: DomainFacetSpec<TColumn, TId, TKind>[];
}

/** One category bucket of a `'categories'` facet: `total` ignores the filter, `filtered` respects
 * it (minus the facet's own column); ref columns group by id and carry the referenced row's
 * display name in `display`. */
export interface DomainFacetCategory {
  /** The category's value as stored. */
  value: any;
  /** The value's display name, where the column has one (a ref's name column). */
  display?: string;
  /** Rows in the category, the filter ignored. */
  total: number;
  /** Rows in the category under the filter. */
  filtered: number;
}
/** Result of a `'categories'` facet; `hasMore` set when the category cap was hit. */
export interface DomainFacetCategoriesResult {
  /** The categories, most frequent first. */
  categories: DomainFacetCategory[];
  /** Whether categories past the cap were left out. */
  hasMore?: boolean;
}
/** Result of a `'histogram'` facet: `buckets` respect the filter, `totalBuckets` ignore it. */
export interface DomainFacetHistogramResult {
  /** The axis start, the filter ignored. */
  min: number | string | null;
  /** The axis end, the filter ignored. */
  max: number | string | null;
  /** Row counts per bucket under the filter. */
  buckets: number[];
  /** Row counts per bucket, the filter ignored. */
  totalBuckets: number[];
  /** Rows without a value. */
  nulls: number;
}
/** Result of a `'minMax'` facet (row predicate only — the stable-axis rule ignores the filter). */
export interface DomainFacetMinMaxResult {
  /** The smallest value. */
  min: number | string | null;
  /** The largest value. */
  max: number | string | null;
}
/** Result of a `'count'` facet: the filtered row count. */
export interface DomainFacetCountResult {
  /** Rows under the filter. */
  count: number;
}
/** Result of a `'plan'` facet: per-column profile for choosing filter-control types. */
export interface DomainFacetPlanResult {
  /** One profile per column asked for: the distinct count (capped) and the range of a numeric or datetime one. */
  columns: {name: string; distinct: number; min?: number | string; max?: number | string}[];
}
/** Maps a facet `kind` to its result type (keys the typed `facets()` response). */
export type DomainFacetResultOf<K extends DomainFacetKind> =
  K extends 'categories' ? DomainFacetCategoriesResult :
  K extends 'histogram' ? DomainFacetHistogramResult :
  K extends 'minMax' ? DomainFacetMinMaxResult :
  K extends 'count' ? DomainFacetCountResult : DomainFacetPlanResult;

/** One operation of `DomainsDataSource.transaction`. */
export interface DomainTransactionOp {
  /** `'restore'` carries `id` alone and undoes a landed soft delete — the Delete
   * grant, not a new permission; the server orders a parent's restore before its
   * child's. */
  op: 'insert' | 'update' | 'delete' | 'restore';
  /** `'<table>'` in the transaction's schema, or `'<schema>.<table>'` — one
   * transaction may span schemas. */
  table: string;
  /** Names this op's new row id; other ops may reference it in values as `'$<ref>'`
   * — forward references are resolved (the server orders the inserts), and deletes
   * run child-first. Escape a literal leading `$` in a value by doubling it:
   * `'$$100'` stores `'$100'`. */
  ref?: string;
  values?: object;
  id?: string;
  /** Optimistic-concurrency guard for update ops on a platform-stored table
   * (`support.concurrency === 'version'`). */
  expectedVersion?: number;
  /** Old-value guard for update ops on an external table (`support.concurrency === 'expected'`):
   * the columns' values as last read, every one AND-ed to the key; a mismatch refuses the whole
   * transaction with a {@link DomainVersionConflictError} carrying `expected` / `current` and the
   * `opIndex`. Never together with {@link expectedVersion}; a literal leading `$` is doubled as in
   * `values`. */
  expected?: {[column: string]: unknown};
  /** Insert ops only: `'error'` fails the whole transaction on a business-key
   * conflict. Without it the insert merges into the existing row and reports
   * `{created: false, status: 'duplicate', existingId}`. */
  onDuplicate?: 'error';
}

/** Options for domain table `batch`. */
export interface DomainBatchOptions {
  /** `'insert'` (default) or `'upsert'` (merge by the table's business key). */
  mode?: 'insert' | 'upsert';
  /** Abort the whole batch on any row error (default true); false applies good rows
   * and reports bad ones per row. */
  allOrNothing?: boolean;
  /** Report business-key duplicates as errors instead of skipping them. */
  errorOnDuplicate?: boolean;
  /** Reject with a {@link DomainValidationError} carrying the report when the batch is aborted. */
  throwOnError?: boolean;
  /** Payload format for `Uint8Array` data: `'d42'` (default; `DataFrame.toByteArray()` output)
   * or `'parquet'` (converted client-side via the Arrow package — fails with a clear error when
   * Arrow is not installed). Ignored for other payloads — the format is inferred:
   * DataFrame → sent as d42, string → `'csv'`, object[] → `'json'`. */
  format?: 'csv' | 'd42' | 'parquet' | 'json';
  /** Judge the payload and write NOTHING: the server runs the whole commit path — coercions,
   * per-row validation, required and auto-number checks, intra-batch and live business-key
   * duplicates, FK existence, immutable columns, and the merge itself — then rolls the
   * transaction back, so the verdicts are the ones a real commit would produce. Resolves to a
   * {@link DomainBatchValidation} instead of a {@link DomainBatchReport}; `allOrNothing` is
   * forced false (a preview judges every row). */
  validateOnly?: boolean;
}

/** Per-row verdict of a `batch({validateOnly: true})` preview. `predicted` is what the row WOULD
 * do; there is no `id` and no `status` — the id of a rolled-back insert does not exist. */
export interface DomainBatchValidationRow {
  index: number;
  predicted: 'insert' | 'update' | 'skip' | 'error';
  /** The row this payload row matched by business key — set for `'update'` and `'skip'` only. */
  existingId?: string;
  errors?: DomainColumnError[];
}

/** What `batch({validateOnly: true})` answers: the counts a commit of this payload would report
 * and the per-row verdicts (capped server-side at 1000, ordered errors → skips → updates →
 * inserts). Its own shape, so a predicted insert can never be read as a completed write.
 *
 * A prediction against the table AS IT IS NOW: a concurrent write between the preview and the
 * commit can change any verdict. Nothing is written, audited, counted or notified by a
 * validation — not a row, not an audit entry, not an auto-number, not a watcher. */
export interface DomainBatchValidation {
  validateOnly: true;
  rowCount: number;
  willInsert: number;
  willUpdate: number;
  willSkip: number;
  errorCount: number;
  rows: DomainBatchValidationRow[];
}

/** Batch upload report of domain table `batch`. */
export interface DomainBatchReport {
  inserted: number;
  updated: number;
  /** Rows an upsert landed on a storage that cannot tell an insert from an update (an external
   * binding): counted here, not under `inserted`/`updated`. A platform table never reports it. */
  merged?: number;
  skipped: number;
  errorCount: number;
  /** Per-row outcomes; capped server-side. */
  rows: DomainBatchRowResult[];
  /** Set when the batch failed but a per-row report is available (e.g. an allOrNothing abort). */
  error?: string;
}

/** The table's change token ({@link DomainTableClient.version}): `seq` moves by one per write
 * TRANSACTION that touched the table's rows (not per row), `at` is when it last moved — null
 * until the first write. */
export interface DomainTableVersion {
  seq: number;
  at: string | null;
}

/** Insert payload for domain table `insert`: row values plus an optional
 * idempotency key (for tables that declare `"idempotency": true`). */
export type DomainRowInsert<TRow> = Partial<TRow> & {idempotencyKey?: string};

/**
 * One link of an expanded many-to-many relation: the target row's id and its display
 * name (name column → dash-joined business key → id, exactly as the row renders
 * everywhere else). A relation expand (`expand: ['labels']`) returns these arrays under
 * the relation's own name, capped at 100 and ordered by display name; the array is `[]`,
 * never null, for an owner with no visible links.
 *
 * `queryDf` flattens the same array into TWO flat string columns instead — `labels`, the
 * display names joined by `', '` and tagged so the grid draws chips, and its companion
 * `'~labels.id'`, the ids in the same order. **The ids column is the source of truth**
 * (a display name containing the separator is sanitized in the flat column, never in the
 * JSON shape); an owner with no links has no value in both — which in a string column of
 * a DataFrame IS the empty string (`col.isNone(i)` is true and `col.get(i)` is `''`;
 * there is no separate null slot), so test emptiness, not `null`.
 *
 * Filtering goes through the same name: `'labels.name'`, `'labels.id'` (a list compiles
 * to ANY-of), and chains that continue from the target (`'labels.group_id.name'`) — one
 * hop each, per-hop EXISTS semantics, so `labels.name != 'bug'` means "has a label that
 * is not bug", not "has no bug label". The exceptions are both on the `'<relation>.id'`
 * leaf, which asks about the link SET rather than about one link — it is what the relation
 * facet's checkboxes emit: `= null` selects owners with NO visible link and `!= null`
 * those with any (the "(no value)" bucket), and `!= [id, ...]` EXCLUDES the owners linked
 * to any of them (the uncheck gesture), rather than "has some other link".
 * Values are bound server-side as everywhere.
 * The smart-filter string form takes the same paths (`'labels.name = "bug"'`); only a list
 * of ids needs the condition tree.
 *
 * Writing is the inverse: `insert` and `update` take the relation as a list of target
 * ids ({@link DomainTableClient.insert}).
 */
export interface DomainRelationLink {
  id: string;
  name: string;
}

/** A saved filter preset of a domain table — a small shareable entity carrying the filter
 * panel's state maps verbatim (see `DomainSavedFiltersClient`). */
export interface DomainSavedFilterInfo {
  id: string;
  name: string;
  friendlyName: string;
  /** Per-column filter states (the shape `DG.FilterGroup` saves), stored verbatim. */
  states: {[column: string]: any};
  /** The preset's author (preserved on in-place updates). */
  author?: any;
}

/** Options of `DomainsDataSource.table`. Datetime columns resolve from the domain
 * registry by default — these overrides exist for legacy generated clients and for
 * callers that must avoid the registry. */
export interface DomainTableClientOptions {
  /** OVERRIDE: datetime columns to materialize as dayjs on JSON reads, instead of the
   * registry-resolved set. Dotted `'<fk>.<col>'` entries cover master-expand fields. */
  datetimeColumns?: string[];
  /** OVERRIDE: datetime columns of `'details:'` child rows, keyed by the result field
   * (the child-table name), instead of the registry-resolved set. */
  detailDatetimeColumns?: {[detailField: string]: string[]};
}

/** Resolved datetime columns of a table: its own plus those of its `'details:'`
 * child tables, keyed by the result field (see `DomainTableClient`). */
export interface DomainDatetimeColumns {
  own: string[];
  details: {[detailField: string]: string[]};
}

/** Typed transaction-op values: each column also accepts a `'$<ref>'` back-reference. */
export type DomainTxValues<T> = {[K in keyof T]: T[K] | `$${string}`};

/** Result type of one transaction op, keyed on its `op` discriminant — powers the
 * mapped-tuple `transaction()` signatures (per-op result types, no positional casts). */
export type DomainOpResultFor<TOp> =
  TOp extends {op: 'insert'} ? DomainInsertResult :
  TOp extends {op: 'update'} ? DomainUpdateResult :
  TOp extends {op: 'delete'} ? DomainDeleteResult :
  TOp extends {op: 'restore'} ? DomainRestoreResult : DomainOpResult;

/** Runs [action], retrying on DomainVersionConflictError (for transaction-based
 * read-modify-write flows — put the fresh read INSIDE [action]); rethrows anything else
 * and the final conflict. `maxRetries` counts retries AFTER the initial attempt
 * (default 5 retries = up to 6 attempts). Retries run immediately, without backoff —
 * conflicts resolve by re-reading, not waiting; add your own delay for high-contention
 * hot rows. */
export async function retryOnVersionConflict<T>(
    action: () => Promise<T>, options?: {maxRetries?: number}): Promise<T> {
  const maxRetries = options?.maxRetries ?? 5;
  for (let attempt = 0; ; attempt++) {
    try {
      return await action();
    } catch (e) {
      if (attempt >= maxRetries || !(e instanceof DomainVersionConflictError))
        throw e;
    }
  }
}

/** Minimal structural client the builder executes against (avoids a domains.ts → dapi.ts cycle). */
export interface IDomainQueryExecutor<TRow> {
  query(spec: any): Promise<TRow[]>;
  queryDf(spec: any): Promise<any>;       // DG.DataFrame
  count(filter?: any, options?: {search?: string}): Promise<number>;
  /** Table address, needed by {@link DomainQueryBuilder.toQuery} (`DomainTableClient` carries it). */
  readonly schema?: string;
  readonly table?: string;
}

/** Bound condition node: `cond('name', '=', "O'Brien")`. The value travels in the condition
 * tree and is bound server-side — NEVER interpolated into a filter string, so any string
 * value is expressible (apostrophes included, which the smart-filter grammar cannot quote). */
export function cond<TColumn extends string = string>(
  property: DomainColumnRef<TColumn>, operator: DomainConditionOperator,
  value?: DomainFilterValue): DomainCondition<TColumn> {
  return value === undefined ? {property, operator} : {property, operator, value};
}

/** AND-combined condition tree: `and(a, b, c)` → `[a, 'and', b, 'and', c]`. */
export function and<TColumn extends string = string>(
  ...nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[]): DomainConditionTree<TColumn> {
  return _joinNodes(nodes, 'and');
}

/** OR-combined condition tree: `or(a, b, c)` → `[a, 'or', b, 'or', c]`. */
export function or<TColumn extends string = string>(
  ...nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[]): DomainConditionTree<TColumn> {
  return _joinNodes(nodes, 'or');
}

function _joinNodes<TColumn extends string>(
  nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[],
  connector: 'and' | 'or'): DomainConditionTree<TColumn> {
  const tree: DomainConditionTree<TColumn> = [];
  for (const n of nodes) {
    if (tree.length > 0)
      tree.push(connector);
    tree.push(n as any);
  }
  return tree;
}

/** Thenable query builder (returned by no-arg `query()`): build with
 * `.where/.orderBy/.select/.expand/.top/.skip`, then `await` it (rows), or finish with
 * `.df()/.first()/.count()/.exists()`. Prefer the condition forms of `where` (and the
 * `cond`/`and`/`or` helpers) over template-built filter strings — condition values are
 * bound server-side, so any string value is safe (apostrophes included). Without `.top()`
 * the server's default limit (100) applies; page larger sets with `.top()/.skip()`.
 * Immutable-ish: `expand()`/`select()` return a re-typed builder; other methods mutate
 * and return this. Invariants: the builder is PromiseLike only — there is no
 * `.catch()`/`.finally()`, use `try { await b } catch`; EVERY `await` (or terminal)
 * re-executes the query — two awaits are two round trips, cache the rows instead. */
export class DomainQueryBuilder<TRow, TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}, TResult = TRow, TDf = any>
    implements PromiseLike<TResult[]> {
  private _conds: DomainConditionTree<TColumn> = [];
  private _rawFilter?: string;
  private _orders: string[] = [];
  private _columns?: string[];
  private _expand: string[] = [];
  private _limit?: number;
  private _offset?: number;
  private _withAccess = false;
  private _search?: string;

  constructor(private readonly client: IDomainQueryExecutor<TRow>) {}

  /** Equality map (AND-combined), a single typed condition (3-arg or node), or a raw
   * smart-filter string escape hatch; multiple where() calls AND-combine. A raw string
   * cannot be combined with conditions (the string parses server-side) — use
   * `cond()/and()/or()` to express everything as one tree instead. */
  where(values: {[K in TColumn]?: DomainFilterValue}): this;
  where(property: DomainColumnRef<TColumn>, operator: DomainConditionOperator,
        value?: DomainFilterValue): this;
  where(filter: DomainFilter<TColumn>): this;
  where(a: any, operator?: DomainConditionOperator, value?: DomainFilterValue): this {
    if (typeof a === 'string' && operator !== undefined)
      this._addCond(cond(a, operator, value) as DomainCondition<TColumn>);
    else if (typeof a === 'string') {
      if (this._rawFilter !== undefined || this._conds.length > 0)
        throw new Error('cannot combine a string filter with conditions — express everything as one tree via cond()/and()/or()');
      this._rawFilter = a;
    }
    else if (Array.isArray(a) || (a != null && 'property' in a && 'operator' in a))
      this._addCond(a);
    else if (a != null)
      for (const k of Object.keys(a))
        this._addCond({property: k as any, operator: '=', value: a[k]});
    return this;
  }

  private _addCond(node: DomainCondition<TColumn> | DomainConditionTree<TColumn>): void {
    if (this._rawFilter !== undefined)
      throw new Error('cannot combine a string filter with conditions — express everything as one tree via cond()/and()/or()');
    if (this._conds.length > 0)
      this._conds.push('and');
    this._conds.push(node as any);
  }

  /** Appends a sort key (the `'col,!col2'` grammar; [desc] adds the `!`). */
  orderBy(column: DomainColumnRef<TColumn>, desc: boolean = false): this {
    this._orders.push(`${desc ? '!' : ''}${column}`);
    return this;
  }

  /** Narrows the projection to [columns] (system columns always ride along).
   * NB: select() re-types the result from TRow — call it BEFORE expand(): the
   * select-then-expand order composes, while expand-then-select keeps the expanded
   * fields at runtime but drops them from the compile-time type. */
  select<K extends TColumn & keyof TRow & string>(...columns: K[]):
      DomainQueryBuilder<TRow, TColumn, TExpand,
        Pick<TRow, K | Extract<keyof TRow, DomainSystemColumn>>, TDf> {
    this._columns = columns;
    return this as any;
  }

  /** Adds an expand and intersects its fields into the awaited row type
   * (`'details:'` child rows are full child rows — dayjs datetimes included). */
  expand<K extends keyof TExpand & string>(key: K):
      DomainQueryBuilder<TRow, TColumn, TExpand, TResult & TExpand[K], TDf> {
    this._expand.push(key);
    return this as any;
  }

  /** Row cap for this query (server default 100 without it). */
  top(count: number): this {
    this._limit = count;
    return this;
  }

  /** Row offset (pair with {@link top} to page). */
  skip(count: number): this {
    this._offset = count;
    return this;
  }

  /** Adds the per-row {@link DOMAIN_ACCESS_COLUMNS} (see {@link DomainQuerySpec.withAccess}). */
  withAccess(): this {
    this._withAccess = true;
    return this;
  }

  /** Substring search over the table's searchable columns (see {@link DomainQuerySpec.search}). */
  search(text: string): this {
    this._search = text;
    return this;
  }

  private _filter(): DomainFilter<TColumn> | undefined {
    return this._rawFilter !== undefined ? this._rawFilter
      : this._conds.length === 0 ? undefined : this._conds;
  }

  private _spec(): any {
    const spec: any = {};
    const filter = this._filter();
    if (filter !== undefined)
      spec.filter = filter;
    if (this._orders.length > 0)
      spec.sort = this._orders.join(',');
    if (this._columns != null)
      spec.columns = this._columns;
    if (this._expand.length > 0)
      spec.expand = this._expand;
    if (this._limit != null)
      spec.limit = this._limit;
    if (this._offset != null)
      spec.offset = this._offset;
    if (this._withAccess)
      spec.withAccess = true;
    if (this._search != null && this._search !== '')
      spec.search = this._search;
    return spec;
  }

  private _run(): Promise<TResult[]> {
    return this.client.query(this._spec()) as Promise<any>;
  }

  /** The same query as a typed DataFrame (d42; `'details:'` expand is JSON-only and
   * rejected server-side — await the builder for detail arrays instead). */
  df(): Promise<TDf> {
    return this.client.queryDf(this._spec());
  }

  /** First matching row or null. NB: permanently overwrites the builder's limit with 1 —
   * a later `await` of the same builder returns at most one row. */
  async first(): Promise<TResult | null> {
    this._limit = 1;
    const rows = await this._run();
    return rows.length === 0 ? null : rows[0];
  }

  /** Matching-row count under the built filter and search (ignores top/skip/select/expand). */
  count(): Promise<number> {
    return this.client.count(this._filter(), {search: this._search});
  }

  /** Whether at least one row matches. */
  async exists(): Promise<boolean> {
    return (await this.count()) > 0;
  }

  /** This builder's accumulated state as a serializable {@link DomainQuery} (the
   * URL / deep-link / recorded-run form of the same query). Conditions become filter
   * elements: top-level AND conjuncts are emitted one per element (so a URL can bind
   * `filters[0]` alone), anything else travels as one JSON element.
   *
   * The row cap travels, but its DEFAULT does not: without `.top()`, awaiting the
   * builder takes the server default of 100 rows while `toQuery().toSpec()` falls back
   * to {@link DOMAIN_QUERY_ROW_LIMIT} — the same builder, two row counts. Pin one with
   * `.top()` when the two forms must agree.
   *
   * A {@link search} is REFUSED rather than dropped: `DomainQuery` is the `DomainQuery`
   * function's parameters, and the function takes no search — a query that carried one
   * would silently select more rows than the collection it came from. Express the same
   * narrowing as a condition, or keep the search as the UI state it is. */
  toQuery(): DomainQuery {
    if (this.client.schema == null || this.client.table == null)
      throw new Error('the query builder has no table address — construct the DomainQuery explicitly');
    if (this._search != null && this._search !== '')
      throw new Error(`toQuery() cannot carry the search "${this._search}": a DomainQuery has no ` +
        'search parameter. Drop the search, or express it as a filter condition.');
    return new DomainQuery({
      schema: this.client.schema, table: this.client.table,
      filters: this._rawFilter !== undefined ? [this._rawFilter] : _treeToFilterElements(this._conds),
      columns: this._columns, joins: this._expand, orderBy: this._orders,
      limit: this._limit, offset: this._offset,
    });
  }

  then<TR1 = TResult[], TR2 = never>(
    onfulfilled?: ((value: TResult[]) => TR1 | PromiseLike<TR1>) | null,
    onrejected?: ((reason: any) => TR2 | PromiseLike<TR2>) | null): Promise<TR1 | TR2> {
    return this._run().then(onfulfilled, onrejected);
  }
}

/** The `DomainQuery` function's client-side row ceiling (the server clamps to it):
 * the limit {@link DomainQuery.toSpec} writes when the query declares none, so a spec
 * never silently inherits the server's default of 100. */
export const DOMAIN_QUERY_ROW_LIMIT = 10000000;

/** List-typed `DomainQuery` parameters, in the order {@link DomainQuery.toUrlParams} emits them. */
const _DOMAIN_QUERY_LISTS = ['columns', 'filters', 'joins', 'aggregations', 'groupBy', 'orderBy'] as const;
type _DomainQueryList = (typeof _DOMAIN_QUERY_LISTS)[number];

function _isListParam(name: string): name is _DomainQueryList {
  return (_DOMAIN_QUERY_LISTS as readonly string[]).includes(name);
}

/** Copy of a list parameter; an empty or absent list normalizes to undefined, so
 * every DomainQuery has exactly ONE representation of "no elements". */
function _copyList(list?: string[] | null): string[] | undefined {
  return list == null || list.length === 0 ? undefined : list.slice();
}

/** Condition tree → `filters` elements: top-level AND conjuncts become one element each
 * (a URL can then bind `filters[0]` alone); anything with an 'or' at the top travels as
 * a single `'['`-prefixed JSON sub-group, which the function splices verbatim. */
function _treeToFilterElements(tree: DomainConditionTree): string[] | undefined {
  if (tree == null || tree.length === 0)
    return undefined;
  if (tree.some((n) => n === 'or'))
    return [JSON.stringify(tree)];
  return tree.filter((n) => typeof n !== 'string').map((n) => JSON.stringify(n));
}

/** One `filters` element as a condition node, or undefined when it is a smart-filter
 * grammar string (which only the function parses — client-side, via the Dart parser). */
function _decodeFilterElement(element: string): any {
  const s = `${element}`.trim();
  if (s.startsWith('{')) {
    let node: any;
    try {
      node = JSON.parse(s);
    } catch (_) {
      throw new Error(`Cannot parse filter "${element}": invalid JSON`);
    }
    if (node == null || typeof node !== 'object')
      throw new Error(`Cannot parse filter "${element}": expected a condition node`);
    return node;
  }
  // A '['-element that decodes as JSON is a condition sub-group; one that does not
  // is grammar ('[bracketed column] = ...').
  if (s.startsWith('[')) {
    let decoded: any;
    try {
      decoded = JSON.parse(s);
    } catch (_) { /* falls through to the grammar */ }
    if (Array.isArray(decoded)) {
      // Same minimal shape check the function applies (domain_query_func.dart): members
      // are conditions, nested groups, or 'and'/'or' — a bare-string member would splice
      // server-side nonsense.
      if (decoded.length === 0 || !decoded.every((e: any) =>
        (e != null && typeof e === 'object') || e === 'and' || e === 'or'))
        throw new Error(`Cannot parse filter "${element}": expected a condition sub-group`);
      return decoded;
    }
  }
  return undefined;
}

/** `limit`/`offset` from a URL: a negative value is rejected rather than passed on —
 * the server clamps it to 0, which would silently return nothing. */
function _parseUrlInt(key: string, value: string): number {
  const s = `${value}`.trim();
  if (!/^\d+$/.test(s))
    throw new Error(`Malformed URL parameter '${key}': expected a non-negative integer, got '${value}'`);
  return parseInt(s, 10);
}

/**
 * What the user is looking at, as ONE serializable object: the parameters of the
 * platform's `DomainQuery` function (filters, joins, aggregations, groupBy, orderBy,
 * projection, limit, offset). URL routing, deep links, saved filters, "open in Table
 * View" and data-synced dashboards all serialize exactly this — there is no parallel
 * query-state vocabulary.
 *
 * ```ts
 * const q = new DG.DomainQuery({schema: 'grit', table: 'issue',
 *   filters: ['status = "open"'], orderBy: ['!created_on'], limit: 100});
 * const df = await q.run();                                  // recorded run
 * const url = new URLSearchParams(q.toUrlParams()).toString(); // filters[0]=...&orderBy[0]=...
 * ```
 *
 * **UI-only state stays out of it.** Search text, view mode, the current entity ride
 * separate reserved URL parameters (`view=`, `entity=`) that {@link fromUrlParams}
 * ignores; if an app needs to carry them together, wrap a DomainQuery in an envelope
 * object — never widen this class.
 */
export class DomainQuery {
  schema: string;
  table: string;
  /** Projection; omit for all viewable columns (select mode only). */
  columns?: string[];
  /** Filter elements, AND-combined: smart-filter grammar strings (`'status = "open"'`),
   * or the per-element JSON escape hatch (a `'{'`-prefixed condition node, a
   * `'['`-prefixed condition sub-group). Values inside a JSON node are bound
   * server-side; a grammar string cannot quote apostrophes, so prefer nodes for
   * arbitrary user values. */
  filters?: string[];
  /** Master-FK expands (`'<fk_column>'`); `'details:'` child arrays are rejected. */
  joins?: string[];
  /** Measures (`'count'`, `'avg(amount) as avg_amount'`) — non-empty means aggregate mode. */
  aggregations?: string[];
  /** Grouping columns — non-empty means aggregate mode. */
  groupBy?: string[];
  /** Sort elements; a `'!'` prefix means descending. */
  orderBy?: string[];
  /** Row cap; {@link toSpec} falls back to {@link DOMAIN_QUERY_ROW_LIMIT} when unset. */
  limit?: number;
  /** Row offset (select mode only). */
  offset?: number;

  /** Empty and absent lists are the same thing: they normalize to undefined, so
   * {@link toParams} / {@link toUrlParams} round-trip to an identical object. */
  constructor(params: DomainQueryParams) {
    if (params == null || params.schema == null || params.table == null)
      throw new Error("DomainQuery needs a 'schema' and a 'table'");
    this.schema = params.schema;
    this.table = params.table;
    for (const name of _DOMAIN_QUERY_LISTS)
      this[name] = _copyList(params[name]);
    if (params.limit != null)
      this.limit = params.limit;
    if (params.offset != null)
      this.offset = params.offset;
  }

  /** The query behind a view's current state: `DG.DomainQuery.fromParams(domainView.query)`. */
  static fromParams(params: DomainQueryParams): DomainQuery {
    return new DomainQuery(params);
  }

  /** The state a {@link DomainQueryBuilder} accumulated (see its `toQuery`). */
  static fromBuilder(builder: DomainQueryBuilder<any>): DomainQuery {
    return builder.toQuery();
  }

  /** The `DomainQuery` function's parameter values — the shape {@link DomainView.query}
   * reports and {@link run} passes to the function. */
  toParams(): DomainQueryParams {
    const params: DomainQueryParams = {schema: this.schema, table: this.table};
    for (const name of _DOMAIN_QUERY_LISTS) {
      const list = _copyList(this[name]);
      if (list !== undefined)
        params[name] = list;
    }
    if (this.limit != null)
      params.limit = this.limit;
    if (this.offset != null)
      params.offset = this.offset;
    return params;
  }

  /** Aggregate mode: `aggregations` or `groupBy` carries elements. `columns` and
   * `offset` are then REJECTED, not ignored — {@link toSpec} throws and {@link run}
   * gets a 400 from the function. */
  get isAggregate(): boolean {
    return (this.aggregations?.length ?? 0) > 0 || (this.groupBy?.length ?? 0) > 0;
  }

  /** Runs the query through the platform's `DomainQuery` function — the reproducible
   * path: the resulting DataFrame carries a creation script, so it refreshes from the
   * Source pane, survives a saved project as a data-synced dashboard, and takes URL
   * parameters (`filters[0]`, ...). Being a normal function run, its result also goes
   * through the platform's default handling: the frame is added to the workspace and
   * OPENED IN A TABLE VIEW, which becomes the current view. Per-caller row and column
   * security applies on every run.
   *
   * The creation script is recorded only when the user has data history enabled
   * (`grok.shell.settings.dataHistory`, on by default) — the query itself runs either
   * way, but the frame comes back without `.script`/`.history` tags when it is off.
   *
   * Failures arrive as the function's own error (raw Dart message text), NOT as the
   * typed {@link DomainError} family: the function boundary carries no error envelope.
   * Where typed failures matter, run {@link toSpec} through `queryDf` instead.
   *
   * For a silent read (no history, no workspace entry, no view) pass {@link toSpec} to
   * `grok.dapi.domains.table('<schema>.<table>').queryDf(...)` instead. */
  async run(): Promise<DataFrame> {
    const func = DG.Func.byName('DomainQuery');
    if (func == null)
      throw new Error("the 'DomainQuery' function is not registered on this client");
    const call = await func.prepare(this.toParams()).call(false, undefined, {processed: false});
    return call.getOutputParamValue();
  }

  /** The same query as a REST spec for `queryDf`/`query` — the same rows as {@link run},
   * without the recording (the limit is always explicit, so the result never depends on
   * the server default).
   *
   * Throws for an aggregate query (use {@link run}, or `aggregateDf` with a
   * {@link DomainAggregateSpec}), and for several smart-filter grammar elements at
   * once: those are parsed by the function itself, and a bare string inside a REST
   * condition tree means a connector, not a filter. */
  toSpec(): DomainQuerySpec {
    if (this.isAggregate)
      throw new Error('an aggregate DomainQuery has no queryDf spec — use run(), or ' +
        'grok.dapi.domains.table(...).aggregateDf() with a DomainAggregateSpec');
    const spec: DomainQuerySpec = {limit: this.limit ?? DOMAIN_QUERY_ROW_LIMIT, offset: this.offset ?? 0};
    const filter = this._toFilter();
    if (filter !== undefined)
      spec.filter = filter;
    if (this.orderBy != null)
      spec.sort = this.orderBy.join(',');
    if (this.joins != null)
      spec.expand = this.joins;
    if (this.columns != null)
      spec.columns = this.columns;
    return spec;
  }

  private _toFilter(): DomainFilter | undefined {
    const elements = (this.filters ?? []).filter((e) => `${e}`.trim() !== '');
    if (elements.length === 0)
      return undefined;
    const nodes: DomainConditionTree = [];
    for (const element of elements) {
      const node = _decodeFilterElement(element);
      if (node === undefined) {
        if (elements.length > 1)
          throw new Error(`toSpec() cannot combine ${elements.length} filter elements into one REST ` +
            `filter: "${element}" is a smart-filter string, and only the DomainQuery function parses ` +
            'those. Use run(), or express the filters as condition nodes (cond()/and()/or()).');
        return `${element}`.trim();
      }
      if (nodes.length > 0)
        nodes.push('and');
      nodes.push(node);
    }
    return nodes;
  }

  /** URL parameters for a deep link, binding one list element per key
   * (`filters[0]=status %3D "open"`) — the platform's list-element binding scheme, so a
   * recorded run's `filters[0]` can be substituted from the URL. The table address
   * (`schema`/`table`) is NOT emitted: it addresses the view, and {@link fromUrlParams}
   * takes it back explicitly. Values are raw — encode them (`URLSearchParams`). */
  toUrlParams(): {[key: string]: string} {
    const params: {[key: string]: string} = {};
    for (const name of _DOMAIN_QUERY_LISTS) {
      const list = this[name];
      if (list != null)
        for (let i = 0; i < list.length; i++)
          params[`${name}[${i}]`] = list[i];
    }
    if (this.limit != null)
      params['limit'] = `${this.limit}`;
    if (this.offset != null)
      params['offset'] = `${this.offset}`;
    return params;
  }

  /** Rebuilds the query from {@link toUrlParams} output (lossless round trip). Keys that
   * are not query parameters — the reserved `view=` / `entity=` UI state among them —
   * are ignored; a malformed element index (`filters[x]`) or a non-integer
   * `limit`/`offset` throws instead of degrading to NaN. Element indices are read in
   * ascending order and gaps are closed. */
  static fromUrlParams(schema: string, table: string, params: {[key: string]: string}): DomainQuery {
    const query = new DomainQuery({schema, table});
    const lists: {[name: string]: {index: number, value: string}[]} = {};
    for (const key of Object.keys(params ?? {})) {
      const value = params[key];
      const bracket = key.indexOf('[');
      const name = bracket < 0 ? key : key.substring(0, bracket);
      if (bracket < 0) {
        if (name === 'limit')
          query.limit = _parseUrlInt(key, value);
        else if (name === 'offset')
          query.offset = _parseUrlInt(key, value);
        else if (_isListParam(name))
          throw new Error(`Malformed URL parameter '${key}': '${name}' is a list — ` +
            `address its elements as '${name}[0]'`);
        continue;
      }
      if (!_isListParam(name))
        continue;
      const index = /^\[(\d+)]$/.exec(key.substring(bracket));
      if (index == null)
        throw new Error(`Malformed URL parameter '${key}': expected '${name}[<index>]'`);
      if (lists[name] == null)
        lists[name] = [];
      lists[name].push({index: parseInt(index[1], 10), value: value});
    }
    for (const name of Object.keys(lists))
      query[name as _DomainQueryList] = lists[name].sort((a, b) => a.index - b.index).map((e) => e.value);
    return query;
  }
}
