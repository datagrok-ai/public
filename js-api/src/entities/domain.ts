/**
 * Domain schema entities: registered, entity-mapped PostgreSQL schemas that hold
 * application data (plates, studies, compounds...). See {@link grok.dapi.domains}.
 * @module entities/domain
 */

import {toJs} from '../wrappers';
import {IDartApi} from '../api/grok_api.g';
import {Entity} from './entity';

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
