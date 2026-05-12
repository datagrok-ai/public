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

// ---------- property-panel side: status, duration, metadata ----------

export function buildExecutionMeta(state: NodeExecState): HTMLElement {
  const container = ui.div([], 'funcflow-value-inspector');
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

// ---------- docked-panel side: just the rich previews ----------

/** True if this state has at least one output worth rendering as a preview
 *  (DataFrame, Column with a sample, or graphics). Used by the docked panel
 *  to decide whether to open at all — primitive-only states stay in the
 *  property panel. */
export function hasRenderablePreview(state: NodeExecState): boolean {
  if (!state.outputs) return false;
  for (const summary of Object.values(state.outputs)) {
    if (summary.type === 'dataframe' && summary.clone) return true;
    if (summary.type === 'column' && Array.isArray(summary.sample) && summary.sample.length > 0) return true;
    if (summary.type === 'graphics' && typeof summary.value === 'string') return true;
  }
  return false;
}

export function buildValuePreviews(state: NodeExecState): HTMLElement {
  const container = ui.div([], 'funcflow-value-previews');
  if (!state.outputs) return container;
  for (const [name, summary] of Object.entries(state.outputs)) {
    const preview = buildPreview(name, summary);
    if (preview) container.appendChild(preview);
  }
  return container;
}

/** A single rich preview for one output value. Returns null for primitives /
 *  null / object — those don't merit dedicating space in the docked panel.
 *  Exported because the per-port "View output" popup uses it directly. */
export function buildPreview(name: string, summary: ValueSummary): HTMLElement | null {
  switch (summary.type) {
  case 'dataframe': {
    if (!summary.clone) return null;
    const wrap = ui.div([], 'funcflow-preview-block');
    try {
      summary.clone.meta.detectSemanticTypes();
      const grid = DG.Viewer.grid(summary.clone as DG.DataFrame);
      grid.root.style.cssText = 'width:100%;min-height:350px;';
      wrap.appendChild(grid.root);
    } catch { /* grid render failed — show nothing rather than a placeholder */ }
    return wrap;
  }
  case 'column': {
    if (!summary.sample || summary.sample.length === 0) return null;
    const wrap = ui.div([], 'funcflow-preview-block');
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
    const wrap = ui.div([], 'funcflow-preview-block');
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
  default:
    return null;
  }
}
