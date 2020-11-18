/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { ChordViewer } from './chord-viewer.js';
import '../css/chord-viewer.css';

export const _package = new DG.Package();

//name: Chord
//description: Creates a chord diagram
//tags: viewer
//output: viewer result
export function chord() {
    return new ChordViewer();
}
