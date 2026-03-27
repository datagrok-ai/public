import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MpxProject} from './mpx-types';
import {mwxWorksheetToDataFrame} from '../mwx/mwx-to-dataframe';


/** Builds a preview for a parsed MPX project with tabs for each worksheet. */
export function buildMpxView(project: MpxProject): HTMLElement {
  // Project info header
  const headerItems: HTMLElement[] = [];
  if (project.creator)
    headerItems.push(ui.divText(`Creator: ${project.creator}`));
  if (project.comments)
    headerItems.push(ui.divText(project.comments));
  const header = ui.div(headerItems, {style: {padding: '8px', borderBottom: '1px solid var(--grey-2)'}});

  // Tab control: one tab per worksheet with lazy grid initialization
  const tabControl = ui.tabControl();
  for (const ws of project.worksheets) {
    tabControl.addPane(ws.name, () => {
      const df = mwxWorksheetToDataFrame(ws, 'mpx');
      const grid = DG.Viewer.grid(df);
      grid.root.style.width = '100%';
      grid.root.style.height = '100%';
      return grid.root;
    });
  }

  if (project.history.length > 0) {
    tabControl.addPane('Session History', () => {
      const d = ui.div([], {style: {whiteSpace: 'pre-wrap', fontSize: '12px', overflow: 'auto', padding: '8px'}});
      d.textContent = project.history.join('\n');
      return d;
    });
  }
  tabControl.root.style.width = '100%';
  tabControl.root.style.flex = '1';

  const container = ui.divV([header, tabControl.root]);
  container.style.width = '100%';
  container.style.height = '100%';
  container.style.display = 'flex';
  container.style.flexDirection = 'column';
  return container;
}


/** Opens all worksheets from an MPX project as table views. */
export function openMpxProject(project: MpxProject): void {
  for (const ws of project.worksheets) {
    const df = mwxWorksheetToDataFrame(ws, 'mpx');
    grok.shell.addTableView(df);
  }
}
