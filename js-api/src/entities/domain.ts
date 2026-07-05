/**
 * Domain schema entities: registered, entity-mapped PostgreSQL schemas that hold
 * application data (plates, studies, compounds...). See {@link grok.dapi.domains}.
 * @module entities/domain
 */

import {toJs} from '../wrappers';
import {IDartApi} from '../api/grok_api.g';
import {Entity} from './entity';
import dayjs from 'dayjs';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** A registered domain schema declared by a plugin manifest (databases/<schema>/schema.json). */
export class DomainSchema extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Physical PostgreSQL schema name. */
  get pgSchema(): string { return api.grok_DomainSchema_Get_PgSchema(this.dart); }

  /** 'package' | 'user' | 'orphaned'. */
  get managedBy(): string { return api.grok_DomainSchema_Get_ManagedBy(this.dart); }

  /** Manifest version applied at last deploy. */
  get version(): string { return api.grok_DomainSchema_Get_Version(this.dart); }

  /** Tables registered in this schema (populate with `.include('tables')`). */
  get tables(): DomainTable[] { return toJs(api.grok_DomainSchema_Get_Tables(this.dart)); }
}


/** A securable table inside a {@link DomainSchema}; its rows are typed by a
 * per-table entity type named '<schema>.<table>'. */
export class DomainTable extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get schema(): DomainSchema { return toJs(api.grok_DomainTable_Get_Schema(this.dart)); }

  /** Row security mode: 'table' | 'master' | 'row'. */
  get securityMode(): string { return api.grok_DomainTable_Get_SecurityMode(this.dart); }

  /** Natural key columns: power dedup-on-insert, upsert matching, and entity handles. */
  get businessKey(): string[] { return toJs(api.grok_DomainTable_Get_BusinessKey(this.dart)); }

  /** Whether writes leave an in-transaction audit trail. */
  get audit(): boolean { return api.grok_DomainTable_Get_Audit(this.dart); }
}


/** A single row of a {@link DomainTable}. Its semantic type is the row's entity
 * type name, `'<schema>.<table>'` — the string an {@link ObjectHandler} keys on to
 * override rendering per table (see Grit's `grit.issue` handler). */
export class DomainRow {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Domain schema name (also the physical PostgreSQL schema). */
  get schemaName(): string { return api.grok_DomainRow_Get_SchemaName(this.dart); }

  /** Table name within the schema. */
  get tableName(): string { return api.grok_DomainRow_Get_TableName(this.dart); }

  /** Row entity type and semantic type: `'<schema>.<table>'`. */
  get typeName(): string { return api.grok_DomainRow_Get_TypeName(this.dart); }

  /** Display/addressing identity: business-key values joined by `'-'`, or the row id. */
  get semValue(): string { return api.grok_DomainRow_Get_SemValue(this.dart); }

  /** Raw row values keyed by wire column name (declared columns, jsonb keys, system fields). */
  get values(): {[key: string]: any} { return toJs(api.grok_DomainRow_Get_Values(this.dart)); }

  /** Row id (GUID) — the storage and security identity. */
  get id(): string { return api.grok_DomainRow_Get_Id(this.dart); }

  /** Optimistic-concurrency version. */
  get version(): number { return api.grok_DomainRow_Get_Version(this.dart); }

  /** Id of the user who last authored the row. */
  get authorId(): string { return api.grok_DomainRow_Get_AuthorId(this.dart); }

  get createdOn(): dayjs.Dayjs | null {
    const d = api.grok_DomainRow_Get_CreatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  get updatedOn(): dayjs.Dayjs | null {
    const d = api.grok_DomainRow_Get_UpdatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  toString(): string { return this.semValue; }
}
