import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { LAYOUT_STORAGE } from './constants';

export type RevvityLayout = {
    tags: {[key: string]: [string, string][]},
    look: {[key: string]: string}
}

export function saveRevvityLayout(grid: DG.Grid, libName: string, entityType: string) {
    const look = grid.getOptions().look;
    if (look.isGrid == undefined)
        look.isGrid = true;
    const tags: {[key: string]: [string, string][]} = {};
    for (const colName of grid.dataFrame.columns.names()) {
        tags[colName] = Object.entries(grid.dataFrame.col(colName)?.tags);
    }
    const revvityLayout: RevvityLayout = { look, tags };
    grok.userSettings.add(LAYOUT_STORAGE, `${libName}|${entityType}`, JSON.stringify(revvityLayout));
}

export async function applyRevvityLayout(grid: DG.Grid, layoutKey: string) {
    const savedLayout = grok.userSettings.getValue(LAYOUT_STORAGE, layoutKey);
    if (savedLayout) {
        try {
            const layout: RevvityLayout = JSON.parse(savedLayout);
            for (const colName of Object.keys(layout.tags)) {
                layout.tags[colName].forEach((tag) => grid.dataFrame.col(colName)?.setTag(tag[0], tag[1]))
            }
            grid.setOptions(layout.look);
        } catch (e) {
            grok.shell.warning(`Failed to restore saved layout: ${e}`);
        }
    }
}