/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import { SunburstViewer } from './sunburst-viewer';

export const _package = new DG.Package();

//name: Sunburst
//description: Sunburst chart
//tags: viewer
//output: viewer result
export function sunburst() {
    return new SunburstViewer();
}
