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

    const dfInp = ui.input.table('Input matrix', {value: [this.solution]});
    dfInp.onChanged(async () => {
      try {
        const df = dfInp.value;        
        dfInp.root.hidden = true;
        div.classList.remove('ui-form');
        const cols = df.columns;
        const timeItems = cols.names().slice(INIT_COL_IDX, -1);
        const timesCount = timeItems.length;
        const names = cols.byIndex(INPUT_NAMES_IDX).toList() as string[];

        const initCol = cols.byIndex(INIT_COL_IDX);
        let arr = initCol.getRawData() as Float32Array;

        const valuesInps = new Map<string, DG.InputBase[]>();

        let currentTimeIdx = 0;

        const updateInputs = () => valuesInps.forEach((inps, time) => {
          inps.forEach((inp) => inp.root.hidden = time !== timeItems[currentTimeIdx]);
        });

        const updateBtns = () => {
          forwardBtn.disabled = currentTimeIdx === timesCount - 1;
          backBtn.disabled = currentTimeIdx === 0;
        };

        const updateCtrls = () => {
          updateInputs();
          updateBtns();
          timeInput.value = timeItems[currentTimeIdx];
        };

        const timeInput = ui.input.string('', {value: timeItems[currentTimeIdx]});
        timeInput.input.style.alignContent = 'center';

        const forwardBtn = ui.button('>', () => {
          if (currentTimeIdx < timesCount - 1) {
            ++currentTimeIdx;
            updateCtrls();
          }
        }, 'Next stage');        

        const backBtn = ui.button('<', () => {
          if (currentTimeIdx > 0) {
            --currentTimeIdx;
            updateCtrls();
          }
        }, 'Previous stage');

        div.append(ui.h2('Time'));
        div.append(ui.divH([backBtn, timeInput.root, forwardBtn]));
        updateBtns();

        const inputsForm = ui.div([ui.h2('Inputs')]);
        timeItems.forEach((time) => {
          arr = cols.byName(time).getRawData() as Float32Array;          
          const inps = names.map((name, idx) => ui.input.float(
            name,
            {
              value: arr[idx],
              onValueChanged: async () => {
                df.set(timeItems[currentTimeIdx], idx, inps[idx].value);
                this.view.dataFrame = await solve(df);                
              },
            }
          ));

          inps.forEach((inp) => inp.root.hidden = time !== cols.byIndex(INIT_COL_IDX).name);

          valuesInps.set(time, inps);

          inputsForm.append(...inps.map((inp) => inp.root));
        });
        inputsForm.classList.add('ui-form');
        inputsForm.style.overflowY = 'scroll';
        div.append(inputsForm);
        
        this.view.dataFrame = await solve(df);
      } catch (e) {
        grok.shell.error(`Failed to solve: ${e instanceof Error ? e.message : 'platform issue'}`);
      }
    });

    div.append(dfInp.root);
    div.classList.add('ui-form');
    div.style.overflowY = 'scroll';
    this.view.dockManager.dock(div, DG.DOCK_TYPE.LEFT, null, undefined, 0.25);    
  }
}