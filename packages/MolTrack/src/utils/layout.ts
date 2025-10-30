import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { LAYOUT_STORAGE, Scope } from './constants';

export type MolTrackLayout = {
    tags: {[key: string]: [string, string][]},
    look: {[key: string]: string}
}

export function saveMolTrackLayout(grid: DG.Grid, scope: Scope) {
  const look = grid.getOptions().look;
  if (look.isGrid == undefined)
    look.isGrid = true;
  const tags: {[key: string]: [string, string][]} = {};
  for (const colName of grid.dataFrame.columns.names())
    tags[colName] = Object.entries(grid.dataFrame.col(colName)?.tags);

  const moltrackLayout: MolTrackLayout = { look, tags };
  //TODO: uncomment user settings when size of storage is increased
  //grok.userSettings.add(LAYOUT_STORAGE, scope, JSON.stringify(moltrackLayout));
  window.localStorage.setItem(`${LAYOUT_STORAGE}|${scope}`, JSON.stringify(moltrackLayout));
}

export async function applyMolTrackLayout(grid: DG.Grid, scope: Scope) {
  //TODO: uncomment user settings when size of storage is increased
  //const savedLayout = grok.userSettings.getValue(LAYOUT_STORAGE, scope);
  const savedLayout = window.localStorage.getItem(`${LAYOUT_STORAGE}|${scope}`);
  if (savedLayout) {
    try {
      const layout: MolTrackLayout = JSON.parse(savedLayout);
      for (const colName of Object.keys(layout.tags))
        layout.tags[colName].forEach((tag) => grid.dataFrame.col(colName)?.setTag(tag[0], tag[1]));

      grid.setOptions(layout.look);
    } catch (e) {
      grok.shell.warning(`Failed to restore saved layout: ${e}`);
    }
  }
}
