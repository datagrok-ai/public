import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solve} from './solver';
import {STAGE, INPUT_NAMES_IDX, INIT_COL_IDX} from './constants';

export class BioRxSim {
  private solution: DG.DataFrame;
  private view: DG.TableView;
  private viewer: DG.Viewer;

  constructor() {
    this.solution = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [0]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'f', [0]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, '_Stage', ['1']),
    ]);
    this.view = grok.shell.addTableView(this.solution);
    this.viewer = DG.Viewer.lineChart(this.solution, {
      segmentColumnName: STAGE,
      multiAxis: true,
    });
    this.view.dockManager.dock(this.viewer, DG.DOCK_TYPE.FILL);

    this.view.name = 'BioRx simulation'
  }

  public async run() {
    const div = ui.div([]);

    const dfInp = ui.input.table('Inputs');
    dfInp.onChanged(async () => {
      try {
        const df = dfInp.value;
        const inpCount = df.rowCount;        
        dfInp.root.hidden = true;
        const cols = df.columns;
        const timeItems = cols.names().slice(INIT_COL_IDX, -1);
        const names = cols.byIndex(INPUT_NAMES_IDX).toList() as string[];

        const initCol = cols.byIndex(INIT_COL_IDX);
        let arr = initCol.getRawData() as Float32Array;

        const timeInp = ui.input.choice<string>('Time', {
          items: timeItems,
          value: initCol.name,
          nullable: false,
          onValueChanged: () => {
            valuesInps.forEach((inps, time) => {
              inps.forEach((inp) => inp.root.hidden = time !== timeInp.value);
            })
          },
        });

        div.append(timeInp.root);

        const valuesInps = new Map<string, DG.InputBase[]>();

        timeItems.forEach((time) => {
          arr = cols.byName(time).getRawData() as Float32Array;          
          const inps = names.map((name, idx) => ui.input.float(
            name,
            {
              value: arr[idx],
              onValueChanged: async () => {
                df.set(timeInp.value, idx, inps[idx].value);
                this.view.dataFrame = await solve(df);                
              },
            }
          ));

          inps.forEach((inp) => inp.root.hidden = time !== cols.byIndex(INIT_COL_IDX).name);

          valuesInps.set(time, inps);

          div.append(...inps.map((inp) => inp.root));
        });                
        
        this.view.dataFrame = await solve(df);
      } catch (e) {
        grok.shell.error(`Failed to solve: ${e instanceof Error ? e.message : 'platform issue'}`);
      }
    });

    div.append(dfInp.root);
    div.classList.add('ui-form');
    this.view.dockManager.dock(div, DG.DOCK_TYPE.LEFT, null, undefined, 0.25);    
  }
}