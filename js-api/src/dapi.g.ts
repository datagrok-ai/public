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
