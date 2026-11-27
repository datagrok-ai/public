/** Configured Input-node values — the run supplies them to the prepared call directly (no parameter
 *  dialog). One editor builder serves the node body and the panel, so the two can never drift. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {FlowNode, hostsInlineSketcher, hostsHelmEditor, editorBoxSize,
  EDITOR_BOX_SIZE_PROP, INLINE_SKETCHER_SIZE, INLINE_EDITOR_WIDTH, WiredEditor,
} from '../rete/scheme';
import {isChoicesRefString, literalChoiceItems} from './choice-refs';
import {loadChoicesRefItems} from '../panel/choice-input-processor';
import {setTid} from './test-ids';
import {effectiveFuncInputs} from './func-input-overrides';

/** Types whose configured value must NOT leak into the `//input:` header default —
 *  a table name / file path / JSON blob is not a valid script default literal. */
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

/** Resolve the configured value into what the prepared call expects for the parameter;
 *  `ok: false` + `reason` when missing or unresolvable. */
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
  case 'string_list': case 'column_list': {
    // An editor's value (`{editor: ...}`) is kept as the list itself: its items may hold commas.
    const items = Array.isArray(raw) ? raw.map((x) => String(x)) : s.split(',').map((x) => x.trim()).filter(Boolean);
    if (items.length > 0) return {ok: true, value: items};
    return nullable ? {ok: true, value: null} : noValue;
  }
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

/** Why this node blocks a silent (dialog-less) run, or null; non-input nodes never block. */
export function inputBlockReason(node: FlowNode): string | null {
  if (node.dgNodeType !== 'input') return null;
  const r = resolveInputValue(node);
  return r.ok ? null : `${node.label} "${String(node.properties['paramName'] ?? '')}": ${r.reason}`;
}

export interface InputValueEditor {
  root: HTMLElement;
  /** Re-read the stored value into the DG input — programmatic, never reported as a user edit. */
  sync: () => void;
  /** The underlying DG input (absent on the sketcher variant, which has no InputBase). */
  input?: DG.InputBase;
  /** The embedded chem sketcher — present only on a Sketcher Input's node-body editor. */
  sketcher?: DG.chem.Sketcher;
  /** Releases what outlives the element (a hosted input's subscription); called before it is thrown away. */
  detach?: () => void;
  /** The data it was built from has changed since (a run replaced its tables): rebuild it. */
  stale?: () => boolean;
}

const SCALAR_PROP_TYPES: Record<string, string> = {
  string: DG.TYPE.STRING, int: DG.TYPE.INT, double: DG.TYPE.FLOAT,
  bool: DG.TYPE.BOOL, datetime: DG.TYPE.DATE_TIME, file: DG.TYPE.FILE,
};

/** The `DG.Property` the value editor is built from — carries the node's qualifiers;
 *  null for types edited as plain text or unsupported. */
export function inputValueProperty(node: FlowNode): DG.Property | null {
  const type = SCALAR_PROP_TYPES[node.dgOutputType ?? ''];
  if (!type) return null;
  const options: Record<string, unknown> = {name: 'Value', type};
  const choicesRaw = String(node.properties['choices'] ?? '').trim();
  if (isChoicesRefString(choicesRaw)) {
    // A reference (func call / query) — the single item IS the reference string;
    // the editor builder swaps in the resolved list asynchronously.
    options['choices'] = [choicesRaw];
  } else {
    const choices = literalChoiceItems(choicesRaw);
    if (choices.length > 0) options['choices'] = choices;
  }
  const min = parseFloat(String(node.properties['min'] ?? ''));
  const max = parseFloat(String(node.properties['max'] ?? ''));
  if (!isNaN(min)) options['min'] = min;
  if (!isNaN(max)) options['max'] = max;
  if (node.properties['nullable'] === true) options['nullable'] = true;
  // A semType routes the editor to the registered (type, semType) valueEditor — a Molecule string is sketched.
  const semType = String(node.properties['semType'] ?? '').trim();
  if (semType) options['semType'] = semType;
  return DG.Property.fromOptions(options as never);
}

/** Build the DG value editor for an input node (null when the type has no inline value).
 *  `onUserChange` fires only on a REAL user edit — never on initialization or `sync()`.
 *  `opts.host: 'node'` builds the node-body variant: a Sketcher Input gets a real
 *  inplace sketcher there, while the panel keeps the compact click-to-sketch editor. */
// eslint-disable-next-line complexity
export function buildInputValueEditor(node: FlowNode, onUserChange: () => void,
  opts?: {host?: 'node' | 'panel'}): InputValueEditor | null {
  const type = node.dgOutputType ?? 'dynamic';
  if (type === 'blob') return null;
  const wired = node.dgNodeType === 'input' ? node.editorBridge?.wiredEditor(node.id) : null;
  if (wired) return buildWiredEditor(node, wired, onUserChange, opts?.host ?? 'panel');
  if (opts?.host === 'node' && hostsInlineSketcher(node))
    return buildInlineSketcherEditor(node, onUserChange);

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
    // An uploaded table exists only inside `ui.input.table`, so keep the live reference;
    // only the name survives save/reload (resolved back by `resolveInputValue`).
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
        else if (type === 'string' && v != null && typeof v !== 'string')
          // A semType-routed editor (Helm) reports a rich value object — the string is in stringValue.
          node.properties['defaultValue'] = input?.stringValue ?? '';
        else
          node.properties['defaultValue'] = v as string | number | boolean ?? '';
      }),
    });
    // Choices stored as a REFERENCE resolve asynchronously into the real item list
    // (the editor starts as a one-item combo holding the reference string). Under
    // the `syncing` guard — swapping items must never read as a user edit.
    const choicesRaw = String(node.properties['choices'] ?? '').trim();
    if (isChoicesRefString(choicesRaw) && input instanceof DG.ChoiceInput) {
      const choiceInput = input as DG.ChoiceInput<string>;
      ui.setUpdateIndicator(choiceInput.input, true, 'loading');
      void loadChoicesRefItems(choicesRaw, String(node.properties['choicesConnection'] ?? ''))
        .then((items) => {
          if (!items || items.length === 0) return;
          syncing = true;
          try {
            choiceInput.items = items;
            const want = String(node.properties['defaultValue'] ?? '');
            if (want !== '' && items.includes(want)) choiceInput.stringValue = want;
          } catch {
            // editor torn down meanwhile
          } finally {
            syncing = false;
          }
        })
        .catch((e) => console.error('Flow: failed to load choices for the value editor', e))
        .finally(() => ui.setUpdateIndicator(choiceInput.input, false));
    }
    const applyValue = (): void => {
      const v = node.properties['defaultValue'];
      if (type === 'bool') input.value = v === true || String(v) === 'true';
      else if (empty(v)) input.value = null;
      else if (type === 'datetime') input.value = dayjs(String(v));
      else input.stringValue = String(v);
    };
    let bindRetry: number | null = null;
    syncValue = (): void => {
      if (bindRetry != null) {
        clearInterval(bindRetry);
        bindRetry = null;
      }
      syncing = true;
      try {
        applyValue();
        return;
      } catch {
        // A semType-routed editor (Helm) materializes asynchronously — a set on
        // the still-unbound proxy throws, and without a retry a loaded flow's
        // value never reaches the editor: it renders (and edits from) empty.
      } finally {
        syncing = false;
      }
      if (empty(node.properties['defaultValue'])) return; // nothing to show
      let tries = 0;
      bindRetry = window.setInterval(() => {
        syncing = true;
        try {
          applyValue();
          clearInterval(bindRetry!);
          bindRetry = null;
        } catch {
          if (++tries > 150) {
            clearInterval(bindRetry!);
            bindRetry = null;
          }
        } finally {
          syncing = false;
        }
      }, 200);
    };
  }
  else {
    // column / column_list / string_list / map / dynamic edit as plain strings; the names resolve at run time.
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
        input.stringValue = empty(v) ? '' : Array.isArray(v) ? v.join(', ') :
          typeof v === 'object' ? JSON.stringify(v) : String(v);
      }
      catch { /* leave the editor as-is */ }
      finally { syncing = false; }
    };
  }

  // Amber underline while the value is missing — the same amber as the ribbon bolt's blocked badge.
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
  const editor: InputValueEditor = {root: input.root, sync, input};
  if (opts?.host === 'node' && hostsHelmEditor(node)) return wrapInHelmBox(node, editor);
  return editor;
}

/** Compact-preview height on a Sketcher Input's node body. */
const SKETCHER_PREVIEW_HEIGHT = 90;

/** The node-body editor of a Sketcher Input: a compact molecule preview that
 *  EXPANDS INTO a 500×500 inplace chem sketcher right on the node — no dialog.
 *  Sketcher chrome is not built to be zoomed, so the sketcher exists ONLY at
 *  native zoom: expanding snaps the canvas to zoom 1 and pans the sketcher fully
 *  into view (`focusNodeForEditing`), and a light poll folds it back the moment
 *  the user zooms away (or the control unmounts).
 *  Stores SMILES — a molfile is multiline and would corrupt the emitted
 *  `//input: string molecule = "…"` header line. The sketcher echoes programmatic
 *  sets back through `onChanged` (asynchronously, canonicalized), so edits report
 *  only after a real user interaction — the `touched` flag, reset on every sync. */
function buildInlineSketcherEditor(node: FlowNode, onUserChange: () => void): InputValueEditor {
  const compact = ui.div([], 'ff-sketcher-compact');
  setTid(compact, 'sketcher-compact');
  ui.tooltip.bind(compact, 'Click to sketch the molecule');
  const full = ui.div([], 'ff-sketcher-full');
  setTid(full, 'sketcher-full');
  full.style.display = 'none';
  const root = ui.div([compact, full], 'ff-inline-sketcher');
  setTid(root, 'inline-sketcher');

  let sketcher: DG.chem.Sketcher | null = null;
  let zoomWatch: number | null = null;
  let touched = false;
  let last = JSON.stringify(node.properties['defaultValue'] ?? '');
  const value = (): string => String(node.properties['defaultValue'] ?? '');
  const markTouched = (): void => {
    touched = true;
  };
  for (const t of ['pointerdown', 'keydown', 'paste'])
    root.addEventListener(t, markTouched, true);
  const markMissing = (): void => {
    root.classList.toggle('ff-value-missing', !resolveInputValue(node).ok);
  };

  const preview = document.createElement('canvas');
  preview.className = 'ff-sketcher-preview';
  const renderCompact = (): void => {
    const v = value();
    ui.empty(compact);
    if (v.trim() === '') {
      compact.appendChild(ui.divText('Sketch', 'ff-sketcher-empty'));
      return;
    }
    compact.appendChild(preview);
    const w = Math.max(compact.clientWidth || 190, 120);
    void grok.chem.canvasMol(0, 0, w, SKETCHER_PREVIEW_HEIGHT, preview, v).catch(() => {
      // Chem not loaded yet — show the raw value instead of a blank box.
      ui.empty(compact);
      compact.appendChild(ui.divText(v, 'ff-sketcher-empty'));
    });
  };

  const setValueIntoSketcher = (): void => {
    try {
      const v = value();
      // Not setMolecule — its SMARTS sniff is a synchronous Chem call that throws
      // until Chem loads; a stored molecule value is SMILES or a molfile.
      if (DG.chem.isMolBlock(v)) sketcher!.setMolFile(v);
      else sketcher!.setSmiles(v);
    } catch {
      // leave the sketcher as-is
    }
  };

  const collapse = (): void => {
    if (zoomWatch != null) {
      clearInterval(zoomWatch);
      zoomWatch = null;
    }
    full.style.display = 'none';
    compact.style.display = '';
    root.classList.remove('ff-sketcher-expanded');
    renderCompact();
  };

  const ensureSketcher = (): void => {
    if (sketcher) return;
    sketcher = new DG.chem.Sketcher();
    editor.sketcher = sketcher;
    sketcher.subs.push(ui.onSizeChanged(sketcher.root).subscribe(() => sketcher!.resize()));
    sketcher.onChanged.subscribe(() => {
      if (!touched) return;
      // The molfile keeps the drawing (coordinates, wedges, query features); the header gets it as SMILES.
      const smiles = sketcher!.getSmiles();
      node.properties['defaultValue'] = sketcher!.getMolFile();
      markMissing();
      const now = JSON.stringify(smiles);
      if (now === last) return;
      last = now;
      onUserChange();
    });
    const done = ui.button('Done', collapse);
    setTid(done, 'sketcher-done');
    ui.tooltip.bind(done, 'Close the sketcher and keep the molecule');
    full.append(ui.divH([done], 'ff-sketcher-header'), sketcher.root);
  };

  const expand = (): void => {
    ensureSketcher();
    touched = false;
    setValueIntoSketcher();
    // The sketcher only works at native zoom — snap to 1 and bring it into view
    // (the box, plus room for the card's title/rows around it).
    node.editorBridge?.focusNodeForEditing(node.id, INLINE_SKETCHER_SIZE + 20, INLINE_SKETCHER_SIZE + 90);
    compact.style.display = 'none';
    full.style.display = '';
    root.classList.add('ff-sketcher-expanded');
    // Fold back the moment the user zooms away (or the control unmounts).
    zoomWatch = window.setInterval(() => {
      const k = node.editorBridge?.getZoom?.() ?? 1;
      if (!root.isConnected || Math.abs(k - 1) > 0.001) collapse();
    }, 200);
  };
  compact.addEventListener('click', expand);

  const sync = (): void => {
    touched = false;
    if (sketcher && root.classList.contains('ff-sketcher-expanded')) setValueIntoSketcher();
    else renderCompact();
    last = JSON.stringify(node.properties['defaultValue'] ?? '');
    markMissing();
  };
  const editor: InputValueEditor = {root, sync};
  sync();
  return editor;
}

/** An input node feeding a function input that has its own editor (`{editor: ...}`, or a registered one) is edited
 *  with that editor, built in the function node's context, so the card reads as the function's panel does; the
 *  value is the input's own. On the card it sits in a resizable box, as big as the editor until the user resizes it. */
function buildWiredEditor(node: FlowNode, wired: WiredEditor, onUserChange: () => void,
  host: 'node' | 'panel'): InputValueEditor {
  const ed = wired.factory(wired.param, wired.ctx);
  let last = JSON.stringify(node.properties['defaultValue'] ?? '');
  ed.onChanged = (v): void => {
    if (ed.isValid && !ed.isValid()) return;
    const now = JSON.stringify(v ?? '');
    if (now === last) return;
    last = now;
    node.properties['defaultValue'] = v as never;
    onUserChange();
  };
  const sync = (): void => {
    last = JSON.stringify(node.properties['defaultValue'] ?? '');
    ed.setValue(node.properties['defaultValue']);
  };
  sync();
  const root = host === 'node' ? resizableBox(node, ed.element, 'ff-editor-box', 'content') : ed.element;
  if (host === 'node') setTid(root, 'editor-box');
  // The editor is built from the function node's captured tables (or a notice that there are none yet).
  const target = wired.ctx.node;
  const tables = (): unknown[] => (target?.dgFunc ? effectiveFuncInputs(target.dgFunc) : [])
    .filter((p) => String(p.propertyType) === 'dataframe')
    .map((p) => wired.ctx.table?.(p.name)?.dart ?? null);
  const builtFrom = tables();
  return {root, sync, detach: ed.detach, stale: () => tables().some((t, i) => t !== builtFrom[i])};
}

/** An in-card box with the native resize handle; a user resize persists into `EDITOR_BOX_SIZE_PROP` — cosmetic,
 *  no params-changed report. Without a stored size it has the default (the helm box's), or with `'content'` follows
 *  the size of what it holds, from {@link INLINE_EDITOR_WIDTH} wide (measured: an editor laid out in percentages has
 *  no intrinsic size to fit by CSS) until the user resizes it. */
function resizableBox(node: FlowNode, content: HTMLElement, cls: string, fit?: 'content'): HTMLDivElement {
  const box = ui.div([content], cls);
  const setSize = (w: number, h: number): void => {
    box.style.width = `${w}px`;
    box.style.height = `${h}px`;
  };
  /** The size last set to fit the content; null once the user has resized the box (or it was never fitting). */
  let fitted: {w: number; h: number} | null = null;
  if (node.properties[EDITOR_BOX_SIZE_PROP] != null || fit !== 'content') {
    const size = editorBoxSize(node);
    setSize(size.width, size.height);
  } else {
    // laid out at a context panel's width first; the height (and any wider part) is then what the editor needs
    fitted = {w: INLINE_EDITOR_WIDTH, h: 0};
    box.style.width = `${INLINE_EDITOR_WIDTH}px`;
    let pending = false;
    // The box scrolls (overflow: auto), so its scroll size is everything it holds, overflow included. Content that
    // grows with the box (a 100%-wide part beside another) overflows by the same amount at any size: once a growth
    // leaves the same overflow, that side stops following it.
    let overflow = {w: 0, h: 0};
    const chasing = {w: false, h: false};
    const fitContent = (): void => {
      pending = false;
      if (!fitted || !box.isConnected || box.scrollWidth === 0) return;
      const now = {w: box.scrollWidth - box.clientWidth, h: box.scrollHeight - box.clientHeight};
      chasing.w ||= now.w > 0 && now.w === overflow.w;
      chasing.h ||= now.h > 0 && now.h === overflow.h;
      overflow = now;
      const cs = getComputedStyle(box);
      const px = (v: string): number => parseFloat(v) || 0;
      const w = chasing.w ? fitted.w :
        Math.max(fitted.w, Math.ceil(box.scrollWidth + px(cs.borderLeftWidth) + px(cs.borderRightWidth)));
      const h = chasing.h ? fitted.h :
        Math.max(fitted.h, Math.ceil(box.scrollHeight + px(cs.borderTopWidth) + px(cs.borderBottomWidth)));
      if (w === fitted.w && h === fitted.h) return;
      fitted = {w, h};
      setSize(w, h);
    };
    const schedule = (): void => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(fitContent);
    };
    // The editor builds (asynchronously) and lays itself out inside: watch its parts as they come.
    const sizes = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(schedule) : null;
    const watch = (): void => {
      for (const el of [content, ...Array.from(content.children)]) sizes?.observe(el);
      schedule();
    };
    new MutationObserver(watch).observe(content, {childList: true, subtree: true});
    watch();
  }
  if (typeof ResizeObserver !== 'undefined') {
    new ResizeObserver(() => {
      if (!box.style.width || !box.style.height) return;
      const w = Math.round(parseFloat(box.style.width));
      const h = Math.round(parseFloat(box.style.height));
      if (!(w > 0) || !(h > 0)) return;
      if (fitted && w === fitted.w && h === fitted.h) return; // the size it was given to fit its content
      fitted = null; // the user resized it: from now on its size is theirs
      const cur = editorBoxSize(node);
      if (cur.width !== w || cur.height !== h)
        node.properties[EDITOR_BOX_SIZE_PROP] = {width: w, height: h};
    }).observe(box);
  }
  return box;
}

/** The Helm Input's node-body editor sits in a resizable in-card box (native CSS
 *  resize handle) so the HELM canvas can be made bigger. The HELM input pins its
 *  editor host to a fixed 250×250 with inline styles and never tracks a
 *  container on its own — stretch `getInput()` to the box (the documented
 *  consumer pattern; the input's own onSizeChanged then re-fits the canvas),
 *  retried until the asynchronously-materialized editor binds. A user resize
 *  persists into `EDITOR_BOX_SIZE_PROP` — cosmetic, no params-changed report. */
function wrapInHelmBox(node: FlowNode, ed: InputValueEditor): InputValueEditor {
  const box = resizableBox(node, ed.root, 'ff-helm-box');
  setTid(box, 'helm-box');
  const fit = (): boolean => {
    let host: HTMLElement | null;
    try {
      host = ed.input?.input ?? null;
    } catch {
      return false; // the JS input proxy hasn't bound yet
    }
    if (!host) return false;
    for (const [k, v] of [['width', '100%'], ['height', '100%'], ['min-width', '0'], ['min-height', '0']])
      host.style.setProperty(k, v, 'important');
    return true;
  };
  if (!fit()) {
    let tries = 0;
    const timer = window.setInterval(() => {
      if (fit() || ++tries > 150) clearInterval(timer);
    }, 200);
  }
  return {root: box, sync: ed.sync, input: ed.input};
}
