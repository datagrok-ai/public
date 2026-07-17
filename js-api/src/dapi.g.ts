// GENERATED CODE - DO NOT MODIFY BY HAND
// Generated from @ExportApi() annotated server plugins
// See core/server/datlas/docs/code-gen.md

let _root: string = '';
let _token: string = '';

export function init(root: string, token: string) {
  _root = root;
  _token = token;
}

async function _fetch(url: string, options: RequestInit = {}): Promise<any> {
  options.headers = { ...options.headers as any, 'Authorization': _token };
  const res = await fetch(`${_root}${url}`, options);
  return res.json();
}

export namespace dapi2 {
  export namespace domains {
    export async function getSchemas(options?: {page?: number, limit?: number, order?: string, include?: string, text?: string}): Promise<any[]> {
      let url = `/domains/schemas`;
      const params = new URLSearchParams();
      if (options?.page !== undefined) params.set('page', String(options.page));
      if (options?.limit !== undefined) params.set('limit', String(options.limit));
      if (options?.order !== undefined) params.set('order', String(options.order));
      if (options?.include !== undefined) params.set('include', String(options.include));
      if (options?.text !== undefined) params.set('text', String(options.text));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    /** Creates an empty user-managed domain schema; body {name, friendlyName?, description?}. The physical PostgreSQL schema is 'usr_<name>'; tables are added afterwards via the apply endpoint. Requires the CreateDomainSchema global privilege; the creator receives the full permission set on the schema entity. */
    export async function createSchema(body: any): Promise<any> {
      let url = `/domains/schemas`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)});
    }

    /** Full manifest of a registered schema, reconstructed from the registry (feeds the schema editor; doubles as export). */
    export async function getSchemaManifest(schema: string): Promise<any> {
      let url = `/domains/schemas/${schema}/manifest`;
      return _fetch(url);
    }

    /** Applies a partial manifest to a user-managed schema: body {tables?, propertySchemas?, dropTables?, ifVersion?, confirmDestructive?}. Named tables replace their current definition wholesale; untouched tables come from the registry. [dryRun] returns the change plan (with live row counts for anything dropped) without applying; a destructive plan requires confirmDestructive; a stale ifVersion is a 409. Requires Edit on the schema entity. */
    export async function applySchema(schema: string, body: any, options?: {dryRun?: boolean}): Promise<any> {
      let url = `/domains/schemas/${schema}/apply`;
      const params = new URLSearchParams();
      if (options?.dryRun !== undefined) params.set('dryRun', String(options.dryRun));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)});
    }

    /** Deletes (fully purges) a user-managed schema: the PostgreSQL schema with its data and audit partition, registry rows, entity types, promoted row entities, and permissions. Requires Delete on the schema entity. */
    export async function deleteSchema(schema: string): Promise<any> {
      let url = `/domains/schemas/${schema}`;
      return _fetch(url, {method: 'DELETE'});
    }

    /** Direct permission rows on a registry entity (a domain schema, table, or property schema): [{group: {id, friendlyName, personal}, permission}]. Requires Share on the entity. */
    export async function getGrants(entityId: string): Promise<any> {
      let url = `/domains/grants/${entityId}`;
      return _fetch(url);
    }

    /** Idempotently grants body 'permission' (View|Edit|Delete|Share) on a registry entity to body 'group' (a group id — sharing with a user resolves to their personal group on the client, as everywhere). Requires Share on the entity. */
    export async function grantPermission(entityId: string, body: any): Promise<any> {
      let url = `/domains/grants/${entityId}`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)});
    }

    /** Revokes 'permission' from 'group' on a registry entity; omitting 'permission' revokes all four common permissions for that group. Requires Share on the entity. */
    export async function revokePermission(entityId: string, options?: {group?: string, permission?: string}): Promise<any> {
      let url = `/domains/grants/${entityId}`;
      const params = new URLSearchParams();
      if (options?.group !== undefined) params.set('group', String(options.group));
      if (options?.permission !== undefined) params.set('permission', String(options.permission));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url, {method: 'DELETE'});
    }

    /** Output format is selected by the spec's 'format': JSON rows (default) or a d42 binary DataFrame (pattern: TablesRouter.getTableData); unknown formats are rejected by the spec parse (400). */
    export async function queryRows(schema: string, table: string, spec: any): Promise<any> {
      let url = `/domains/${schema}/${table}/query`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(spec)});
    }

    export async function aggregateRows(schema: string, table: string, spec: any): Promise<any> {
      let url = `/domains/${schema}/${table}/aggregate`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(spec)});
    }

    /** Bulk upload (§5.6); the body format is selected by Content-Type: application/json (row-object array, bare or under 'rows'), text/csv, or application/octet-stream (d42 DataFrame). No server-side Parquet — see PLAN.md phase-2 decision 1 (clients convert via the Arrow package). The generated dapi2 client method is a non-functional stub (it sends no body): use grok.dapi.domains.table(...).batch (JS) or raw HTTP with a JSON/CSV/d42 body instead. */
    export async function batchRows(schema: string, table: string): Promise<any> {
      let url = `/domains/${schema}/${table}/batch`;
      return _fetch(url, {method: 'POST'});
    }

    /** Ordered per-op results, shaped like the single-op endpoints; any failure rolls the whole transaction back and the error body carries 'opIndex'. [ops] is the ops list; a {'ops': [...]} wrapper is also accepted over the wire (`@RequestBody() dynamic`, not `List` — see insertRows). */
    export async function runTransaction(schema: string, ops: any): Promise<any> {
      let url = `/domains/${schema}/transaction`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(ops)});
    }

    /** [rows] is a single row map or a list of row maps. The body is `@RequestBody() dynamic`, not `Map` — a checked-mode VM enforces handler annotations on the mirror invoke, and the wire accepts both shapes (pinned by domain_rest_test). */
    export async function insertRows(schema: string, table: string, rows: any): Promise<any> {
      let url = `/domains/${schema}/${table}`;
      return _fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(rows)});
    }

    /** Current user's subscription state for the table — {'watching': bool}. */
    export async function getTableWatch(schema: string, table: string): Promise<any> {
      let url = `/domains/${schema}/${table}/watch`;
      return _fetch(url);
    }

    export async function watchTable(schema: string, table: string): Promise<any> {
      let url = `/domains/${schema}/${table}/watch`;
      return _fetch(url, {method: 'POST'});
    }

    export async function unwatchTable(schema: string, table: string): Promise<any> {
      let url = `/domains/${schema}/${table}/watch`;
      return _fetch(url, {method: 'DELETE'});
    }

    export async function getRowWatch(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}/watch`;
      return _fetch(url);
    }

    /** Row-level watch requires the table's audit trail (`audit: false` → 400). */
    export async function watchRow(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}/watch`;
      return _fetch(url, {method: 'POST'});
    }

    export async function unwatchRow(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}/watch`;
      return _fetch(url, {method: 'DELETE'});
    }

    /** Table-wide audit events, newest first; [limit] clamps to [1, 1000]. */
    export async function tableAudit(schema: string, table: string, options?: {limit?: number}): Promise<any> {
      let url = `/domains/${schema}/${table}/audit`;
      const params = new URLSearchParams();
      if (options?.limit !== undefined) params.set('limit', String(options.limit));
      if (params.toString()) url += '?' + params.toString();
      return _fetch(url);
    }

    export async function getRow(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}`;
      return _fetch(url);
    }

    export async function deleteRow(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}`;
      return _fetch(url, {method: 'DELETE'});
    }

    export async function promoteRow(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}/promote`;
      return _fetch(url, {method: 'POST'});
    }

    export async function rowAudit(schema: string, table: string, id: string): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}/audit`;
      return _fetch(url);
    }

    export async function patchRow(schema: string, table: string, id: string, body: any): Promise<any> {
      let url = `/domains/${schema}/${table}/${id}`;
      return _fetch(url, {method: 'PATCH', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)});
    }

  }

  export namespace chats {
  }

}
