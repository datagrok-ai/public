/** Two views over a node's runtime state, used by different surfaces:
 *
 *  - {@link buildExecutionMeta} — status badge, duration, error/stack trace,
 *    and *metadata* about each output (DataFrame dims + column names, column
 *    name + length, etc). Inline values for primitives. Goes into the
 *    property panel ("Execution" section).
 *
 *  - {@link buildValuePreviews} — only the rich previews: DataFrame grids
 *    with "Add to workspace", column sample tables, graphics images. No
 *    headers, no metadata. Goes into the bottom-docked output panel.
 *
 *  Both consume the same {@link NodeExecState}. Splitting them lets the
 *  docked panel show "just the data" while the property panel shows "what
 *  happened" — the user asked to keep the docked panel uncluttered. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NodeExecState, NodeExecStatus, ValueSummary} from './execution-state';
import {setTid} from '../utils/test-ids';

// ---------- property-panel side: status, duration, metadata ----------

export function buildExecutionMeta(state: NodeExecState): HTMLElement {
  const container = setTid(ui.div([], 'funcflow-value-inspector'), 'value-inspector');
  container.appendChild(buildStatusBadge(state));

  if (state.outputs && Object.keys(state.outputs).length > 0) {
    const header = ui.divText('Outputs');
    header.style.fontWeight = 'bold';
    header.style.marginBottom = '4px';
    container.appendChild(header);
    for (const [name, summary] of Object.entries(state.outputs))
      container.appendChild(buildMetaRow(name, summary));
  }

  if (state.error) {
    const errorHeader = ui.divText('Error');
    errorHeader.style.fontWeight = 'bold';
    errorHeader.style.color = '#F44336';
    errorHeader.style.marginTop = '8px';
    container.appendChild(errorHeader);
    const errorMsg = ui.divText(state.error);
    errorMsg.style.color = '#F44336';
    errorMsg.style.fontSize = '12px';
    errorMsg.style.whiteSpace = 'pre-wrap';
    errorMsg.style.wordBreak = 'break-word';
    container.appendChild(errorMsg);

    if (state.stack) {
      const stackDetails = document.createElement('details');
      const stackSummary = document.createElement('summary');
      stackSummary.textContent = 'Stack trace';
      stackSummary.style.cssText = 'cursor:pointer;font-size:11px;color:#999;';
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
    [NodeExecStatus.idle]: '#78909c',
    [NodeExecStatus.running]: '#1976d2',
    [NodeExecStatus.completed]: '#43a047',
    [NodeExecStatus.errored]: '#e53935',
    [NodeExecStatus.stale]: '#9E9E9E',
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
  dot.style.cssText = `width:10px;height:10px;border-radius:50%;background-color:${colors[state.status] ?? '#888'};display:inline-block;`;
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
      }, 'Add to workspace');
      addBtn.style.cssText = 'cursor:pointer;color:#1976d2;font-size:13px;';
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
    // The live object renders in the docked panel; here just name its kind
    // (was "[object Object]" when it fell through to the generic object case).
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

// ---------- "Add to workspace" — overlay button on a preview block ----------

/** A small "add to workspace" button overlaid in the top-right corner of a
 *  preview block (same corner treatment as the viewer gear — no vertical
 *  footprint). `onClick` opens the value in the platform workspace. When a
 *  gear also occupies the corner (viewer previews), pass `rightPx` to sit it to
 *  the gear's left. */
function addWorkspaceButton(title: string, onClick: () => void, rightPx = 6): HTMLElement {
  const btn = setTid(ui.iconFA('plus-circle', (e: Event) => {
    e.stopPropagation();
    onClick();
  }, title), 'add-to-workspace');
  btn.classList.add('ff-preview-workspace-btn');
  btn.style.right = `${rightPx}px`;
  return btn;
}

/** Add a viewer's data to the workspace: opens a fresh table view over a clone
 *  of the viewer's DataFrame, recreates the viewer on that clone (a *new*
 *  viewer — never the previewed instance) with the same look, and docks it to
 *  the right of the table view. */
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

// ---------- docked-panel side: just the rich previews ----------

/** True if this state has at least one output worth rendering as a preview
 *  (DataFrame, Column with a sample, or graphics). Used by the docked panel
 *  to decide whether to open at all — primitive-only states stay in the
 *  property panel. */
export function hasRenderablePreview(state: NodeExecState): boolean {
  if (!state.outputs) return false;
  for (const summary of Object.values(state.outputs)) {
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
  const entries = Object.entries(state.outputs);
  // When the node's real output is a column, its preview is that column rendered
  // as a one-column DataFrame — so suppress the threaded "<input> (modified)"
  // passthrough table (it's kept in the state for the column picker / inspect,
  // but here it's redundant with, and noisier than, the column itself).
  const hasColumnOutput = entries.some(([, s]) => s.type === 'column');
  for (const [name, summary] of entries) {
    if (hasColumnOutput && summary.type === 'dataframe' && name.endsWith('(modified)')) continue;
    const preview = buildPreview(name, summary, onEditViewer);
    if (preview) container.appendChild(preview);
  }
  return container;
}

/** A single rich preview for one output value. Returns null for primitives /
 *  null / object — those don't merit dedicating space in the docked panel.
 *  Exported because the per-port "View output" popup uses it directly. */
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
      const grid = DG.Viewer.grid(df);
      grid.root.style.cssText = 'width:100%;height: calc(100% - 10px);';
      wrap.appendChild(grid.root);
      wrap.appendChild(addWorkspaceButton('Add table to workspace',
        () => grok.shell.addTableView(df.clone())));
    } catch { /* grid render failed — show nothing rather than a placeholder */ }
    return wrap;
  }
  case 'column': {
    // In-place case: the output column was added by instance to the node's
    // single input table (e.g. Add New Column). Show the *whole* table scrolled
    // to the new column — the column in the context of the data it belongs to —
    // rather than a lone one-column grid. (Captured by __ff_col_summary.)
    if (summary.tableClone && summary.scrollToColumn) {
      const df = summary.tableClone as DG.DataFrame;
      const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
      wrap.style.position = 'relative';
      try {
        df.meta.detectSemanticTypes();
        const grid = DG.Viewer.grid(df);
        grid.root.style.cssText = 'width:100%;height: calc(100% - 10px);';
        wrap.appendChild(grid.root);
        // Scroll to the produced column once the grid has laid out (it isn't
        // attached yet here — defer so the horizontal scroll actually applies).
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
    // Preferred: the instrumented run captured a one-column DataFrame clone
    // built from the output column — render it as a real grid.
    if (summary.clone) {
      const df = summary.clone as DG.DataFrame;
      const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
      wrap.style.position = 'relative';
      try {
        df.meta.detectSemanticTypes();
        const grid = DG.Viewer.grid(df);
        grid.root.style.cssText = 'width:100%;height: calc(100% - 10px);';
        wrap.appendChild(grid.root);
        wrap.appendChild(addWorkspaceButton('Add column to workspace',
          () => grok.shell.addTableView(df.clone())));
        return wrap;
      } catch { /* grid failed — fall back to the text sample below */ }
    }
    if (!summary.sample || summary.sample.length === 0) return null;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    const title = ui.divText(`${name}: ${summary.name ?? ''}`);
    title.style.cssText = 'font-size:12px;color:#444;margin-bottom:4px;';
    wrap.appendChild(title);
    const table = document.createElement('table');
    table.style.cssText = 'font-size:11px;color:#555;border-collapse:collapse;';
    const headerRow = table.insertRow();
    const idxTh = document.createElement('th');
    idxTh.textContent = '#';
    idxTh.style.cssText = 'padding:2px 8px;text-align:right;color:#999;font-weight:normal;';
    headerRow.appendChild(idxTh);
    const valTh = document.createElement('th');
    valTh.textContent = 'Value';
    valTh.style.cssText = 'padding:2px 8px;text-align:left;font-weight:normal;color:#999;';
    headerRow.appendChild(valTh);
    for (let i = 0; i < summary.sample.length; i++) {
      const tr = table.insertRow();
      const idxTd = tr.insertCell();
      idxTd.textContent = String(i);
      idxTd.style.cssText = 'padding:1px 8px;text-align:right;color:#999;';
      const valTd = tr.insertCell();
      valTd.textContent = String(summary.sample[i]);
      valTd.style.cssText = 'padding:1px 8px;border-bottom:1px solid #eee;';
    }
    wrap.appendChild(table);
    return wrap;
  }
  case 'graphics': {
    const imageData = summary.value as string;
    if (typeof imageData !== 'string') return null;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    const img = ui.div([], {style: {
      width: '100%', minHeight: '200px',
      backgroundPosition: 'left', backgroundRepeat: 'no-repeat', backgroundSize: 'contain',
      aspectRatio: '1',
    }});
    if (imageData.startsWith('<svg')) img.innerHTML = imageData;
    else img.style.backgroundImage = `url('data:image/png;base64,${imageData}')`;
    wrap.appendChild(img);
    return wrap;
  }
  case 'widget':
  case 'viewer': {
    // The live DG.Widget / DG.Viewer captured during the run — mount its root
    // directly. Its root isn't attached anywhere else, so this just adopts it.
    const obj = summary.value as {root?: HTMLElement} | undefined;
    if (!obj?.root || !(obj.root instanceof Element)) return null;
    const wrap = setTid(ui.div([], 'funcflow-preview-block'), 'preview-block', name);
    wrap.style.position = 'relative';
    obj.root.style.width = '100%';
    if (!obj.root.style.minHeight) obj.root.style.height = 'calc(100% - 10px)';
    wrap.appendChild(obj.root);
    // A viewer's full settings are editable live: a small gear in the top-right
    // corner (overlaid — no vertical space) hands the viewer to the host
    // (→ grok.shell.o = viewer), which captures changes back onto the node.
    const hasGear = summary.type === 'viewer' && !!onEditViewer;
    if (hasGear) {
      const gear = setTid(ui.iconFA('cog', () => onEditViewer!(obj),
        'Edit viewer settings'), 'viewer-edit');
      gear.classList.add('ff-viewer-edit-gear');
      wrap.appendChild(gear);
    }
    // Every viewer is DataFrame-backed — "Add to workspace" recreates it (a new
    // viewer, same look) on a fresh table view. Sits to the gear's left when
    // both are shown, else in the corner.
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
