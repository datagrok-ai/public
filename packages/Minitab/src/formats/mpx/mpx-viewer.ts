import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MpxProject} from './mpx-types';
import {mwxWorksheetToDataFrame} from '../mwx/mwx-to-dataframe';


/** Builds a preview for a parsed MPX project. */
export function buildMpxView(project: MpxProject, onOpen: () => void): HTMLElement {
  const sections: HTMLElement[] = [];

  // Open button
  const openBtn = ui.bigButton('Open', onOpen);
  sections.push(ui.div([openBtn], {style: {padding: '8px'}}));

  // Project info
  const infoItems: HTMLElement[] = [];
  if (project.creator)
    infoItems.push(ui.divText(`Creator: ${project.creator}`));
  if (project.comments)
    infoItems.push(ui.divText(`Comments: ${project.comments}`));
  infoItems.push(ui.divText(`Worksheets: ${project.worksheets.length}`));
  sections.push(ui.div(infoItems, {style: {padding: '8px', borderBottom: '1px solid var(--grey-2)'}}));

  // Worksheets table
  const wsTable = ui.table(
    project.worksheets,
    (ws) => [ws.name, `${ws.rowCount}`, `${ws.columns.length}`],
    ['Worksheet', 'Rows', 'Columns'],
  );
  wsTable.style.width = '100%';
  sections.push(ui.div([ui.h2('Worksheets'), wsTable], {style: {padding: '8px'}}));

  // Session history (collapsible)
  if (project.history.length > 0) {
    const historyText = project.history.join('\n');
    const pre = document.createElement('pre');
    pre.textContent = historyText;
    pre.style.whiteSpace = 'pre-wrap';
    pre.style.fontSize = '12px';
    pre.style.maxHeight = '200px';
    pre.style.overflow = 'auto';
    pre.style.background = 'var(--grey-1)';
    pre.style.padding = '8px';
    pre.style.borderRadius = '4px';

    const details = document.createElement('details');
    const summary = document.createElement('summary');
    summary.textContent = `Session History (${project.history.length} entries)`;
    summary.style.cursor = 'pointer';
    summary.style.fontWeight = 'bold';
    details.appendChild(summary);
    details.appendChild(pre);
    sections.push(ui.div([details], {style: {padding: '8px'}}));
  }

  const container = ui.divV(sections);
  container.style.overflow = 'auto';
  return container;
}


/** Opens all worksheets from an MPX project as table views. */
export function openMpxProject(project: MpxProject): void {
  for (const ws of project.worksheets) {
    const df = mwxWorksheetToDataFrame(ws);
    grok.shell.addTableView(df);
  }
}
