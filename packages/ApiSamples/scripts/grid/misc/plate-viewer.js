// Use grid as a plate viewer

let wells = grok.data.demo.wells(8 * 12);
let plate = wells
  .groupBy(['row'])
  .pivot('col')
  .min('volume')
  .aggregate();

let grid = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, plate);
let view = grok.shell.newView();
view.append(grid);
