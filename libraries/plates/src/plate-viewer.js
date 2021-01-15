/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export function plateViewer() {
  // Use grid as a plate viewer

  let wells = grok.data.demo.wells(8 * 12);
  let plate = wells
    .groupBy(['row'])
    .pivot('col')
    .min('volume')
    .aggregate();

  let grid = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, plate);
  return grid;
}