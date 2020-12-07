/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { SankeyViewer } from '../sankey/sankey.js';
import { GlobeViewer } from '../globe/globe.grok.js';
/* TODO: move cell renderer tests out of this package */
import { FlagCellRenderer } from './flag-cell-renderer.js';

export const _package = new DG.Package();

//name: Sankey
//description: Creates a sankey viewer
//tags: viewer
//output: viewer result
export function sankey() {
    return new SankeyViewer();
}

//name: Globe
//description: Creates a globe viewer
//tags: viewer
//output: viewer result
export function globe() {
    return new GlobeViewer(this.webRoot + '/globe/');
}

//name: flagCellRenderer
//tags: cellRenderer, cellRenderer-flag
//meta-cell-renderer-sem-type: flag
//output: grid_cell_renderer result
export function flagCellRenderer() { 
    return new FlagCellRenderer();
}
