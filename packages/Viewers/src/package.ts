/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';

import {SankeyViewer} from '../sankey/sankey';
import {GlobeViewer} from '../globe/globe-viewer';
import {WordCloudViewer} from '../word-cloud/word-cloud-viewer';
import {ChordViewer} from '../chord/chord-viewer';
import {TreeViewer} from '../tree/tree-viewer';

import '../css/chord-viewer.css';
import '../css/sankey.css';

// /* TODO: move cell renderer tests out of this package */
import {FlagCellRenderer} from './flag-cell-renderer';


export const _package = new DG.Package();

//name: Tree
//description: Phylogenetic tree visualization
//tags: viewer
//output: viewer result
export function tree() {
  return new TreeViewer();
}

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
  return new GlobeViewer();
}

//name: Chord
//description: Creates a chord diagram
//tags: viewer
//output: viewer result
export function chord() {
  return new ChordViewer();
}

//name: Word Cloud
//description: Creates a word cloud
//tags: viewer
//output: viewer result
export function wordcloud() {
  return new WordCloudViewer();
}

//name: flagCellRenderer
//tags: cellRenderer
//meta-cell-renderer-sem-type: flag
//output: grid_cell_renderer result
export function flagCellRenderer() {
  return new FlagCellRenderer();
}
