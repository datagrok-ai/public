/* eslint-disable max-len */
/** Two views over a node's runtime state: {@link buildExecutionMeta} (property
 *  panel — status/metadata) and {@link buildValuePreviews} (output panel — rich values). */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NodeExecState, NodeExecStatus, ValueSummary} from './execution-state';
import {INLINE_HOSTED_DATA_KEY} from '../rete/scheme';
import {setTid} from '../utils/test-ids';
import * as rxjs from 'rxjs';

/** Host callback fired when the user clicks a cell in a preview grid — a
 *  context signal for the suggestion engine. */
let _previewCellFocusHandler:
  ((cell: {semType: string | null; column: string; value: unknown}) => void) | null = null;

export function setPreviewCellFocusHandler(
  h: ((cell: {semType: string | null; column: string; value: unknown}) => void) | null,
): void {
  _previewCellFocusHandler = h;
}

/** Release the hook only if this owner still holds it — a closing view must not
 *  tear down the handler a newer view has since installed. */
export function releasePreviewCellFocusHandler(
  h: ((cell: {semType: string | null; column: string; value: unknown}) => void) | null,
): void {
  if (h !== null && _previewCellFocusHandler === h) {
    _previewCellFocusHandler = null;
    previewHookSub?.unsubscribe();
    previewHookSub = null;
  }
}

let previewHookSub: rxjs.Subscription | null = null;
function hookPreviewCellFocus(df: DG.DataFrame): void {
  previewHookSub?.unsubscribe();
  try {
    previewHookSub = df.onCurrentCellChanged.subscribe(() => {
      if (!_previewCellFocusHandler) return;
      try {
        const cell = df.currentCell;
        const col = cell?.column;
        if (!col) return;
        _previewCellFocusHandler({
          semType: col.semType ? String(col.semType) : null,
          column: col.name,
          value: cell.value,
        });
      } catch {/* cell read failed mid-update — skip this signal */}
    });
  } catch {/* observable unavailable on odd proxies — no signal, no harm */}
}

export function buildExecutionMeta(state: NodeExecState): HTMLElement {
  const container = setTid(ui.div([], 'funcflow-value-inspector'), 'value-inspector');
  container.appendChild(buildStatusBadge(state));

  // `__pt` entries are dims-only bookkeeping for the on-edge count labels (possibly null).
  const shown = Object.entries(state.outputs ?? {})
    .filter(([name, summary]) => summary != null && !name.endsWith('__pt'));
  if (shown.length > 0) {
    // Outputs kept from before a failure must not read as THIS run's result.
    const fromBefore = state.status === NodeExecStatus.errored || state.status === NodeExecStatus.stale;
    const header = ui.divText(fromBefore ? 'Last successful outputs' : 'Outputs');
    header.style.fontWeight = 'bold';
    header.style.marginBottom = '4px';
    container.appendChild(header);
    const rows = ui.div(shown.map(([name, summary]) => buildMetaRow(name, summary)));
    if (fromBefore) rows.style.opacity = '0.6';
    container.appendChild(rows);
  }

  if (state.error) {
    const errorHeader = ui.divText('Error');
    errorHeader.style.fontWeight = 'bold';
    errorHeader.style.color = 'var(--red-3, #eb6767)';
    errorHeader.style.marginTop = '8px';
    container.appendChild(errorHeader);
    const errorMsg = ui.divText(state.error);
    errorMsg.style.color = 'var(--red-3, #eb6767)';
    errorMsg.style.fontSize = '12px';
    errorMsg.style.whiteSpace = 'pre-wrap';
    errorMsg.style.wordBreak = 'break-word';
    container.appendChild(errorMsg);

    if (state.stack) {
      const stackDetails = document.createElement('details');
      const stackSummary = document.createElement('summary');
      stackSummary.textContent = 'Stack trace';
      stackSummary.style.cssText = 'cursor:pointer;font-size:11px;color:var(--grey-4, #9497a0);';
      stackDetails.appendChild(stackSummary);
      const stackPre = document.createElement('pre');
      stackPre.textContent = state.stack;
      stackPre.style.cssText = 'font-size:10px;max-height:150px;overflow:auto;margin:4px 0;';
      stackDetails.appendChild(stackPre);
      container.appendChild(stackDetails);
    }
  }

  return container;
}

function buildStatusBadge(state: NodeExecState): HTMLElement {
  const colors: Record<string, string> = {
    [NodeExecStatus.idle]: 'var(--grey-4, #9497a0)',
    [NodeExecStatus.running]: 'var(--blue-1, #2083d5)',
    [NodeExecStatus.completed]: 'var(--green-2, #3cb173)',
    [NodeExecStatus.errored]: 'var(--red-3, #eb6767)',
    [NodeExecStatus.stale]: 'var(--grey-3, #b8bac0)',
  };
  const labels: Record<string, string> = {
    [NodeExecStatus.idle]: 'Idle',
    [NodeExecStatus.running]: 'Running...',
    [NodeExecStatus.completed]: 'Completed',
    [NodeExecStatus.errored]: 'Error',
    [NodeExecStatus.stale]: 'Stale',
  };

  const badge = ui.div([], 'funcflow-exec-badge');
  badge.style.cssText = 'display:flex;align-items:center;gap:6px;margin-bottom:8px;';

  const dot = document.createElement('span');
  dot.style.cssText = `width:10px;height:10px;border-radius:50%;background-color:${colors[state.status] ?? 'var(--grey-4, #9497a0)'};display:inline-block;`;
  badge.appendChild(dot);

  let label = labels[state.status] ?? state.status;
  if ((state.status === NodeExecStatus.completed || state.status === NodeExecStatus.errored) &&
      state.startTime && state.endTime)
    label += ` (${state.endTime - state.startTime}ms)`;
  badge.appendChild(ui.divText(label));
  return badge;
}

/** Per-output metadata row for the property panel — text only, no grids. */
function buildMetaRow(name: string, summary: ValueSummary): HTMLElement {
  const row = ui.div([], 'funcflow-value-row');
  row.style.cssText = 'margin-bottom:4px;font-size:12px;';

  switch (summary.type) {
  case 'dataframe': {
    const headerLine = ui.div([], {style: {display: 'flex', alignItems: 'center', gap: '6px'}});
    headerLine.appendChild(ui.divText(`${name}: DataFrame (${summary.rows} rows × ${summary.cols} cols)`));
    if (summary.clone) {
      const addBtn = ui.iconFA('plus-circle', () => {
        grok.shell.addTableView(summary.clone as DG.DataFrame);
      }, 'Open this result as a full table view in the workspace');
      addBtn.style.cssText = 'cursor:pointer;color:var(--blue-1, #2083d5);font-size:13px;';
      headerLine.appendChild(addBtn);
    }
    row.appendChild(headerLine);
    if (summary.colNames) {
      const colList = ui.divText(`Columns: ${summary.colNames.join(', ')}`);
      colList.style.cssText = 'font-size:11px;color:#666;margin-left:8px;';
      row.appendChild(colList);
    }
    break;
  }
  case 'column':
    row.appendChild(ui.divText(`${name}: Column "${summary.name}" (${summary.length} values)`));
    break;
  case 'graphics':
    row.appendChild(ui.divText(`${name}: Graphics`));
    break;
  case 'widget':
  case 'viewer':
    row.appendChild(ui.divText(`${name}: ${summary.type}`));
    break;
  case 'primitive':
    row.appendChild(ui.divText(`${name} = ${JSON.stringify(summary.value)}`));
    break;
  case 'object':
    row.appendChild(ui.divText(`${name}: ${summary.str || '[object]'}`));
    break;
  case 'null':
    row.appendChild(ui.divText(`${name} = null`));
    break;
  }
  return row;
}

/** "Add to workspace" button overlaid in a preview block's top-right corner;
 *  pass `rightPx` to sit it left of a gear. */
function addWorkspaceButton(title: string, onClick: () => void, rightPx = 6): HTMLElement {
  const btn = setTid(ui.iconFA('plus-circle', (e: Event) => {
    e.stopPropagation();
    onClick();
  }, title), 'add-to-workspace');
  btn.classList.add('ff-preview-workspace-btn');
  btn.style.right = `${rightPx}px`;
  return btn;
}

/** Recreate the viewer (never the previewed instance) on a clone of its table
 *  in a fresh table view. */
async function addViewerToWorkspace(viewer: DG.Viewer): Promise<void> {
  try {
    const df = viewer.dataFrame;
    if (!df) {
      grok.shell.warning('This viewer has no table to add.');
      return;
    }
    const tv = grok.shell.addTableView(df.clone());
    const opts = viewer.getOptions();
    const newViewer = await tv.dataFrame.plot.fromType(opts.type, {}) as DG.Viewer;
    try {
      newViewer.setOptions(opts.look ?? {});
    } catch {/* look may reference removed columns — leave defaults */}
    tv.dockManager.dock(newViewer, DG.DOCK_TYPE.RIGHT, null, opts.type, 0.4);
  } catch (e) {
    grok.shell.error(`Could not add the viewer to the workspace: ${e instanceof Error ? e.message : e}`);
  }
}

/** Element rendering a graphics output (inline SVG markup or a base64 PNG) —
 *  shared by the bottom panel and the in-node preview. Graphics is data, not a
 *  live object, so each caller builds its own copy. */
export function graphicsElement(imageData: string): HTMLElement {
  const img = ui.div([], {style: {
    width: '100%', minHeight: '200px',
    backgroundPosition: 'left', backgroundRepeat: 'no-repeat', backgroundSize: 'contain',
    aspectRatio: '1',
  }});
  if (imageData.startsWith('<svg')) img.innerHTML = imageData;
  else img.style.backgroundImage = `url('data:image/png;base64,${imageData}')`;
  return img;
}

/** True if this state has at least one output worth a rich preview —
 *  primitive-only states stay in the property panel. */
export function hasRenderablePreview(state: NodeExecState): boolean {
  if (!state.outputs) return false;
  for (const summary of Object.values(state.outputs)) {
    if (summary == null) continue; // `__pt` dims entries can be null
    if (summary.type === 'dataframe' && summary.clone) return true;
    if (summary.type === 'column' && (summary.clone || (Array.isArray(summary.sample) && summary.sample.length > 0))) return true;
    if (summary.type === 'graphics' && typeof summary.value === 'string') return true;
    if ((summary.type === 'widget' || summary.type === 'viewer') && summary.value?.root instanceof Element) return true;
  }
  return false;
}

export function buildValuePreviews(
  state: NodeExecState, onEditViewer?: (viewer: unknown) => void,
): HTMLElement {
  const container = setTid(ui.div([], 'funcflow-value-previews'), 'value-previews');
  if (!state.outputs) return container;
  // `__pt` entries are dims-only bookkeeping for the on-edge count labels (possibly null).
  const entries = Object.entries(state.outputs)
    .filter(([name, summary]) => summary != null && !name.endsWith('__pt'));
  // A column output already previews as a one-column table — the threaded
  // "<input> (modified)" passthrough table would be redundant with it.
  const hasColumnOutput = entries.some(([, s]) => s.type === 'column');
  const blocks: HTMLElement[] = [];
  for (const [name, summary] of entries) {
    if (hasColumnOutput && summary.type === 'dataframe' && name.endsWith('(modified)')) continue;
    const preview = buildPreview(name, summary, onEditViewer);
    if (preview) blocks.push(preview);
  }
  if (blocks.length === 1)
    container.appendChild(blocks[0]);
  else if (blocks.length > 1) {
    const split = ui.splitH(blocks, null, true);
    split.style.width = '100%';
    split.style.height = '100%';
    container.appendChild(split);
  }
  return container;
}

/** A single rich preview for one output value; null for primitives / null / object. */
export function buildPreview(
  name: string, summary: ValueSummary, onEditViewer?: (viewer: unknown) => void,
): HTMLElement | null {
  switch (summary.type) {
  case 'dataframe': {
    if (!summary.clone) return null;
    const df = summary.clone as DG.DataFrame;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    wrap.style.position = 'relative';
    try {
      df.meta.detectSemanticTypes();
      hookPreviewCellFocus(df);
      const grid = DG.Viewer.grid(df);
      grid.root.style.cssText = 'width:100%;height: calc(100% - 16px);';
      wrap.appendChild(grid.root);
      wrap.appendChild(addWorkspaceButton('Add table to workspace',
        () => grok.shell.addTableView(df.clone())));
    } catch {/* grid render failed — show nothing rather than a placeholder */}
    return wrap;
  }
  case 'column': {
    // In-place case (__ff_col_summary): the column was added to the node's input
    // table — show the whole table scrolled to it, not a lone one-column grid.
    if (summary.tableClone && summary.scrollToColumn) {
      const df = summary.tableClone as DG.DataFrame;
      const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
      wrap.style.position = 'relative';
      try {
        df.meta.detectSemanticTypes();
        hookPreviewCellFocus(df);
        const grid = DG.Viewer.grid(df);
        grid.root.style.cssText = 'width:100%;height: calc(100% - 16px);';
        wrap.appendChild(grid.root);
        // The grid isn't attached yet — defer so the horizontal scroll applies.
        const colName = summary.scrollToColumn as string;
        requestAnimationFrame(() => {
          try {
            grid.scrollToCell(colName, 0);
          } catch {/* column gone / grid detached — leave at default scroll */}
        });
        wrap.appendChild(addWorkspaceButton('Add table to workspace',
          () => grok.shell.addTableView(df.clone())));
        return wrap;
      } catch {/* grid failed — fall back to the one-column grid / sample below */}
    }
    if (summary.clone) {
      const df = summary.clone as DG.DataFrame;
      const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
      wrap.style.position = 'relative';
      try {
        df.meta.detectSemanticTypes();
        hookPreviewCellFocus(df);
        const grid = DG.Viewer.grid(df);
        grid.root.style.cssText = 'width:100%;height: calc(100% - 16px);';
        wrap.appendChild(grid.root);
        wrap.appendChild(addWorkspaceButton('Add column to workspace',
          () => grok.shell.addTableView(df.clone())));
        return wrap;
      } catch {/* grid failed — fall back to the text sample below */}
    }
    if (!summary.sample || summary.sample.length === 0) return null;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    const title = ui.divText(`${name}: ${summary.name ?? ''}`);
    title.style.cssText = 'font-size:12px;color:var(--grey-6, #4d5261);margin-bottom:4px;';
    wrap.appendChild(title);
    const table = document.createElement('table');
    table.style.cssText = 'font-size:11px;color:var(--grey-5, #717581);border-collapse:collapse;';
    const headerRow = table.insertRow();
    const idxTh = document.createElement('th');
    idxTh.textContent = '#';
    idxTh.style.cssText = 'padding:2px 8px;text-align:right;color:var(--grey-4, #9497a0);font-weight:normal;';
    headerRow.appendChild(idxTh);
    const valTh = document.createElement('th');
    valTh.textContent = 'Value';
    valTh.style.cssText = 'padding:2px 8px;text-align:left;font-weight:normal;color:var(--grey-4, #9497a0);';
    headerRow.appendChild(valTh);
    for (let i = 0; i < summary.sample.length; i++) {
      const tr = table.insertRow();
      const idxTd = tr.insertCell();
      idxTd.textContent = String(i);
      idxTd.style.cssText = 'padding:1px 8px;text-align:right;color:var(--grey-4, #9497a0);';
      const valTd = tr.insertCell();
      valTd.textContent = String(summary.sample[i]);
      valTd.style.cssText = 'padding:1px 8px;border-bottom:1px solid var(--grey-2, #dbdcdf);';
    }
    wrap.appendChild(table);
    return wrap;
  }
  case 'graphics': {
    const imageData = summary.value as string;
    if (typeof imageData !== 'string') return null;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    wrap.appendChild(graphicsElement(imageData));
    return wrap;
  }
  case 'widget':
  case 'viewer': {
    // The live object captured during the run — its root is attached nowhere else.
    const obj = summary.value as {root?: HTMLElement} | undefined;
    if (!obj?.root || !(obj.root instanceof Element)) return null;
    // A root the in-node preview hosts must not be stolen — one live DOM
    // element cannot be in two places.
    if (obj.root.dataset?.[INLINE_HOSTED_DATA_KEY] === 'true') {
      const note = setTid(ui.div([], 'funcflow-preview-block'), 'preview-inline-note');
      note.classList.add('ff-preview-inline-note');
      note.textContent = `This ${summary.type} is previewed on its node — ` +
        'turn off the node preview (⊟) to show it here.';
      return note;
    }
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    wrap.style.position = 'relative';
    obj.root.style.width = '100%';
    if (!obj.root.style.minHeight) obj.root.style.height = 'calc(100% - 16px)';
    wrap.appendChild(obj.root);
    const hasGear = summary.type === 'viewer' && !!onEditViewer;
    if (hasGear) {
      const gear = setTid(ui.iconFA('cog', () => onEditViewer!(obj),
        'Edit viewer settings'), 'viewer-edit');
      gear.classList.add('ff-viewer-edit-gear');
      wrap.appendChild(gear);
    }
    if (summary.type === 'viewer') {
      const viewer = summary.value as DG.Viewer;
      wrap.appendChild(addWorkspaceButton('Add viewer to workspace',
        () => {void addViewerToWorkspace(viewer);}, hasGear ? 32 : 6));
    }
    return wrap;
  }
  default:
    return null;
  }
}
