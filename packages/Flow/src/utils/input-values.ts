/** Configured values for Input nodes — "set the value right on the node".
 *
 * An input node still compiles to an `//input:` script header (the flow stays
 * a parameterized script), but when a value is configured the run supplies it
 * to the prepared call directly — no parameter dialog, and autorun is not
 * blocked. The value lives in `properties['defaultValue']` (serialized with
 * the flow; for scalar types it doubles as the header default). Values that
 * cannot be expressed as a string — a picked DataFrame, a FileInfo — keep a
 * runtime-only reference in `FlowNode.transientValue` next to the serialized
 * name, so an uploaded table works now and degrades to a by-name lookup after
 * a reload.
 *
 * One editor builder serves both surfaces (the node body control and the
 * context panel), so they can never drift apart; `sync()` re-reads the store,
 * guarded so programmatic updates never count as user edits. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {FlowNode} from '../rete/scheme';

/** Types whose configured value must NOT leak into the `//input:` header
 *  default — a table name / file path / JSON blob is not a valid script
 *  default literal. Scalar types keep the classic `= <value>` emission. */
export const NON_HEADER_DEFAULT_TYPES = new Set(['dataframe', 'file', 'map', 'blob']);

export interface ResolvedInputValue {
  ok: boolean;
  value?: unknown;
  /** Human sentence fragment for the blocked-autorun tooltip ("no table selected"). */
  reason?: string;
}

function empty(v: unknown): boolean {
  return v === undefined || v === null || String(v) === '';
}

/** Resolve the node's configured value into what the prepared script call
 *  expects for its parameter — a live DataFrame for `dataframe`, dayjs for
 *  `datetime`, name strings for columns, a string array for lists. `ok: false`
 *  means the value is missing or unresolvable; `reason` says why. */
// eslint-disable-next-line complexity
export function resolveInputValue(node: FlowNode): ResolvedInputValue {
  const type = node.dgOutputType ?? 'dynamic';
  const raw = node.properties['defaultValue'];
  const s = raw === undefined || raw === null ? '' : String(raw);
  const nullable = node.properties['nullable'] === true;
  const noValue = {ok: false, reason: 'no value set'};

  switch (type) {
  case 'dataframe': {
    if (node.transientValue instanceof DG.DataFrame) return {ok: true, value: node.transientValue};
    if (s === '') return {ok: false, reason: 'no table selected'};
    const df = grok.shell.tables.find((t) => t.name === s) ??
      grok.shell.tables.find((t) => t.name.toLowerCase() === s.toLowerCase());
    return df ? {ok: true, value: df} : {ok: false, reason: `table "${s}" is not open`};
  }
  case 'int': case 'double': {
    if (s === '') return nullable ? {ok: true, value: null} : noValue;
    const n = Number(s);
    if (isNaN(n)) return {ok: false, reason: `"${s}" is not a number`};
    return {ok: true, value: type === 'int' ? Math.round(n) : n};
  }
  case 'bool':
    return {ok: true, value: raw === true || s === 'true'};
  case 'string':
    return s !== '' || nullable ? {ok: true, value: s} : noValue;
  case 'datetime': {
    if (s === '') return nullable ? {ok: true, value: null} : noValue;
    const d = dayjs(s);
    return d.isValid() ? {ok: true, value: d} : {ok: false, reason: `"${s}" is not a valid date`};
  }
  case 'string_list': case 'column_list':
    return s === '' ? noValue : {ok: true, value: s.split(',').map((x) => x.trim()).filter(Boolean)};
  case 'column': case 'dynamic':
    return s === '' ? noValue : {ok: true, value: s};
  case 'map': {
    if (s === '') return noValue;
    try {
      const parsed: unknown = typeof raw === 'object' ? raw : JSON.parse(s);
      return {ok: true, value: parsed};
    }
    catch {
      return {ok: false, reason: 'the value is not valid JSON'};
    }
  }
  case 'file':
    if (node.transientValue != null) return {ok: true, value: node.transientValue};
    return {ok: false, reason: s !== '' ? `pick the file again ("${s}")` : 'no file picked'};
  default: // blob and anything unknown — dialog-only
    return {ok: false, reason: 'this input type has no inline value'};
  }
}

/** Why this node blocks a silent (dialog-less) run, or null when it doesn't.
 *  Non-input nodes never block. */
export function inputBlockReason(node: FlowNode): string | null {
  if (node.dgNodeType !== 'input') return null;
  const r = resolveInputValue(node);
  return r.ok ? null : `${node.label} "${String(node.properties['paramName'] ?? '')}": ${r.reason}`;
}

// ---- value editor (shared by the node body and the context panel) ----------

export interface InputValueEditor {
  root: HTMLElement;
  /** Re-read the stored value into the DG input (programmatic — never reported
   *  as a user edit). Called when the other surface edited the same node. */
  sync: () => void;
}

const SCALAR_PROP_TYPES: Record<string, string> = {
  string: DG.TYPE.STRING, int: DG.TYPE.INT, double: DG.TYPE.FLOAT,
  bool: DG.TYPE.BOOL, datetime: DG.TYPE.DATE_TIME, file: DG.TYPE.FILE,
};

/** The `DG.Property` the value editor is built from — carries the node's
 *  qualifiers (choices, min/max, nullable) so `ui.input.forProperty` renders
 *  the right editor. Null for types edited as plain text or unsupported. */
export function inputValueProperty(node: FlowNode): DG.Property | null {
  const type = SCALAR_PROP_TYPES[node.dgOutputType ?? ''];
  if (!type) return null;
  const options: Record<string, unknown> = {name: 'Value', type};
  const choices = String(node.properties['choices'] ?? '').split(',').map((x) => x.trim()).filter(Boolean);
  if (choices.length > 0) options['choices'] = choices;
  const min = parseFloat(String(node.properties['min'] ?? ''));
  const max = parseFloat(String(node.properties['max'] ?? ''));
  if (!isNaN(min)) options['min'] = min;
  if (!isNaN(max)) options['max'] = max;
  if (node.properties['nullable'] === true) options['nullable'] = true;
  return DG.Property.fromOptions(options as never);
}

/** Build the DG value editor for an input node, or null when the type has no
 *  inline value (blob). `onUserChange` fires only on a REAL user edit — never
 *  on initialization or on `sync()` — mirroring the property panel's
 *  change-reporter guard, so merely rendering a node is not an edit. */
// eslint-disable-next-line complexity
export function buildInputValueEditor(node: FlowNode, onUserChange: () => void): InputValueEditor | null {
  const type = node.dgOutputType ?? 'dynamic';
  if (type === 'blob') return null;

  let syncing = false;
  let last = JSON.stringify(node.properties['defaultValue'] ?? '');
  const report = (): void => {
    const now = JSON.stringify(node.properties['defaultValue'] ?? '');
    if (now === last) return;
    last = now;
    onUserChange();
  };
  const guard = (write: () => void): void => {
    if (syncing) return;
    write();
    report();
  };

  let input: DG.InputBase;
  let syncValue: () => void;

  if (type === 'dataframe') {
    // `ui.input.table` lists the open tables, tracks add/close live, and has a
    // built-in folder icon that opens a LOCAL file into the list — exactly the
    // workspace-or-upload choice this node needs. An uploaded table exists
    // only here, so keep the live reference; the name alone is what survives
    // a save/reload (resolved back from the workspace by `resolveInputValue`).
    const tableInput = ui.input.table('Value', {
      onValueChanged: (v: DG.DataFrame | null) => guard(() => {
        node.transientValue = v ?? undefined;
        node.properties['defaultValue'] = v?.name ?? '';
      }),
    });
    input = tableInput;
    syncValue = (): void => {
      const name = String(node.properties['defaultValue'] ?? '');
      const live = node.transientValue instanceof DG.DataFrame ? node.transientValue : null;
      const df = live ?? (name === '' ? null : grok.shell.tables.find((t) => t.name === name) ?? null);
      try {
        if (df) tableInput.value = df;
        else tableInput.stringValue = name; // shows the saved name even when the table isn't open
      }
      catch { /* leave the editor as-is */ }
    };
  }
  else if (type === 'file') {
    const prop = inputValueProperty(node)!;
    input = ui.input.forProperty(prop, null, {
      onValueChanged: (v: DG.FileInfo | null) => guard(() => {
        node.transientValue = v ?? undefined;
        node.properties['defaultValue'] = v?.fullPath ?? v?.name ?? '';
      }),
    });
    syncValue = (): void => {
      if (node.transientValue == null) return; // a path alone can't rebuild a FileInfo
      try {
        input.value = node.transientValue;
      }
      catch { /* leave the editor as-is */ }
    };
  }
  else if (SCALAR_PROP_TYPES[type]) {
    const prop = inputValueProperty(node)!;
    input = ui.input.forProperty(prop, null, {
      onValueChanged: (v: unknown) => guard(() => {
        if (type === 'datetime')
          node.properties['defaultValue'] = v == null ? '' : (v as dayjs.Dayjs).toISOString();
        else
          node.properties['defaultValue'] = v as string | number | boolean ?? '';
      }),
    });
    syncValue = (): void => {
      const v = node.properties['defaultValue'];
      syncing = true;
      try {
        if (type === 'bool') input.value = v === true || String(v) === 'true';
        else if (empty(v)) input.value = null;
        else if (type === 'datetime') input.value = dayjs(String(v));
        else input.stringValue = String(v);
      }
      catch { /* leave the editor as-is */ }
      finally { syncing = false; }
    };
  }
  else {
    // column / column_list / string_list / map / dynamic — a plain string
    // (comma-separated names, JSON for map). No live table exists to drive a
    // real column picker here; the names resolve when the flow runs.
    const tips: Record<string, string> = {
      column: 'Column name', column_list: 'Column names, comma-separated',
      string_list: 'Values, comma-separated', map: 'JSON, e.g. {"key": "value"}',
      dynamic: 'Value',
    };
    input = ui.input.string('Value', {
      tooltipText: tips[type] ?? 'Value',
      onValueChanged: (v: string) => guard(() => {
        node.properties['defaultValue'] = String(v ?? '');
      }),
    });
    syncValue = (): void => {
      const v = node.properties['defaultValue'];
      syncing = true;
      try {
        input.stringValue = empty(v) ? '' : (typeof v === 'object' ? JSON.stringify(v) : String(v));
      }
      catch { /* leave the editor as-is */ }
      finally { syncing = false; }
    };
  }

  // "Act here" affordances: placeholder text where the editor supports it, and
  // an amber underline while the value is missing/unresolvable — the same
  // amber as the ribbon bolt's blocked badge, so the two cues read as one.
  const placeholders: Record<string, string> = {
    dataframe: 'Choose a table…', string: 'Type a value…', column: 'Column name…',
    column_list: 'Column names, comma-separated…', string_list: 'Values, comma-separated…',
    map: '{"key": "value"}', dynamic: 'Type a value…', file: 'Pick a file…',
  };
  const editorEl = input.input as HTMLInputElement | undefined;
  if (editorEl && 'placeholder' in editorEl && placeholders[type])
    editorEl.placeholder = placeholders[type];
  const markMissing = (): void => {
    input.root.classList.toggle('ff-value-missing', !resolveInputValue(node).ok);
  };
  input.onChanged.subscribe(() => markMissing());

  const sync = (): void => {
    syncing = true;
    try {
      syncValue();
    }
    finally {
      syncing = false;
      last = JSON.stringify(node.properties['defaultValue'] ?? '');
    }
    markMissing();
  };
  sync();
  ui.tooltip.bind(input.root, 'The value used when the flow runs — with it set, ' +
    'Run and autorun don\'t need to ask. Saved with the flow as the parameter default.');
  return {root: input.root, sync};
}
