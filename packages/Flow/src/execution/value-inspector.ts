/** Builds context panel section for displaying runtime execution values */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NodeExecState, NodeExecStatus, ValueSummary} from './execution-state';

/** Creates the execution results section for the property panel */
export function buildValuePanel(state: NodeExecState): HTMLElement {
  const container = ui.div([], 'funcflow-value-inspector');

  // Status badge
  const statusEl = buildStatusBadge(state);
  container.appendChild(statusEl);

  // Duration
  if (state.startTime && state.endTime) {
    const duration = state.endTime - state.startTime;
    const durationEl = ui.divText(`Duration: ${duration}ms`);
    durationEl.style.fontSize = '11px';
    durationEl.style.color = '#666';
    durationEl.style.marginBottom = '8px';
    container.appendChild(durationEl);
  }

  // Outputs
  if (state.outputs && Object.keys(state.outputs).length > 0) {
    const header = ui.divText('Outputs');
    header.style.fontWeight = 'bold';
    header.style.marginBottom = '4px';
    container.appendChild(header);

    for (const [name, summary] of Object.entries(state.outputs))
      container.appendChild(buildValueRow(name, summary));
  }

  // Error
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
      stackSummary.style.cursor = 'pointer';
      stackSummary.style.fontSize = '11px';
      stackSummary.style.color = '#999';
      stackDetails.appendChild(stackSummary);

      const stackPre = document.createElement('pre');
      stackPre.textContent = state.stack;
      stackPre.style.fontSize = '10px';
      stackPre.style.maxHeight = '150px';
      stackPre.style.overflow = 'auto';
      stackPre.style.margin = '4px 0';
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
  badge.style.display = 'flex';
  badge.style.alignItems = 'center';
  badge.style.gap = '6px';
  badge.style.marginBottom = '8px';

  const dot = document.createElement('span');
  dot.style.width = '10px';
  dot.style.height = '10px';
  dot.style.borderRadius = '50%';
  dot.style.backgroundColor = colors[state.status] || '#888';
  dot.style.display = 'inline-block';
  badge.appendChild(dot);

  let label = labels[state.status] || state.status;
  if (state.status === NodeExecStatus.completed && state.startTime && state.endTime)
    label += ` (${state.endTime - state.startTime}ms)`;
  badge.appendChild(ui.divText(label));

  return badge;
}

function buildValueRow(name: string, summary: ValueSummary): HTMLElement {
  const row = ui.div([], 'funcflow-value-row');
  row.style.marginBottom = '4px';
  row.style.fontSize = '12px';

  switch (summary.type) {
  case 'dataframe': {
    // Header line with size info and + button to add to workspace
    const headerLine = ui.div([], {style: {display: 'flex', alignItems: 'center', gap: '6px'}});
    headerLine.appendChild(ui.divText(`${name}: DataFrame (${summary.rows} rows \u00d7 ${summary.cols} cols)`));

    if (summary.clone) {
      const addBtn = ui.iconFA('plus-circle', () => {
        grok.shell.addTableView(summary.clone as DG.DataFrame);
      }, 'Add to workspace');
      addBtn.style.cursor = 'pointer';
      addBtn.style.color = '#1976d2';
      addBtn.style.fontSize = '13px';
      headerLine.appendChild(addBtn);
    }
    row.appendChild(headerLine);

    if (summary.colNames) {
      const colList = ui.divText(`Columns: ${summary.colNames.join(', ')}`);
      colList.style.fontSize = '11px';
      colList.style.color = '#666';
      colList.style.marginLeft = '8px';
      row.appendChild(colList);
    }

    // Full table preview when clone is available
    if (summary.clone) {
      try {
        summary.clone.meta.detectSemanticTypes();
        const grid = DG.Viewer.grid(summary.clone as DG.DataFrame);
        grid.root.style.width = '100%';
        grid.root.style.minHeight = '350px';
        grid.root.style.marginTop = '4px';
        row.appendChild(grid.root);
      } catch { /* grid preview failed, show text only */ }
    }
    break;
  }
  case 'column': {
    row.appendChild(ui.divText(
      `${name}: Column "${summary.name}" (${summary.length} values)`,
    ));
    if (summary.sample && summary.sample.length > 0) {
      const table = document.createElement('table');
      table.style.cssText = 'font-size:11px; color:#555; margin:4px 0 0 8px; border-collapse:collapse;';
      const headerRow = table.insertRow();
      const idxTh = document.createElement('th');
      idxTh.textContent = '#';
      idxTh.style.cssText = 'padding:2px 8px; text-align:right; color:#999; font-weight:normal;';
      headerRow.appendChild(idxTh);
      const valTh = document.createElement('th');
      valTh.textContent = 'Value';
      valTh.style.cssText = 'padding:2px 8px; text-align:left; font-weight:normal; color:#999;';
      headerRow.appendChild(valTh);
      for (let i = 0; i < summary.sample.length; i++) {
        const tr = table.insertRow();
        const idxTd = tr.insertCell();
        idxTd.textContent = String(i);
        idxTd.style.cssText = 'padding:1px 8px; text-align:right; color:#999;';
        const valTd = tr.insertCell();
        valTd.textContent = String(summary.sample[i]);
        valTd.style.cssText = 'padding:1px 8px; border-bottom:1px solid #000000;';
      }
      row.appendChild(table);
    }
    break;
  }
  case 'graphics': {
    row.appendChild(ui.divText(`${name}: Graphics`));
    const imageData = summary.value as string;
    const preview = ui.div([], {style: {
      width: '100%', minHeight: '200px', marginTop: '4px',
      backgroundPosition: 'left', backgroundRepeat: 'no-repeat', backgroundSize: 'contain',
      aspectRatio: '1',
    }});
    if (imageData.startsWith('<svg'))
      preview.innerHTML = imageData;
    else
      preview.style.backgroundImage = `url('data:image/png;base64,${imageData}')`;
    row.appendChild(preview);
    break;
  }
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
