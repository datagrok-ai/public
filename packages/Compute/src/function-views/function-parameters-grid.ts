import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//description: A spreadsheet that lets you interactively edit parameters and evaluate functions
//tags: functionAnalysis
//input: func f
//output: view result
export function _functionParametersGrid(f: DG.Func): DG.View {

  let view = DG.View.create();
  let gridHost = ui.div();
  let functionName = ui.stringInput('Function', f.name);

  view.root.appendChild(functionName.root);
  view.root.appendChild(gridHost);

  function init() {
    view.name = f.name;
    ui.empty(gridHost);

    let table = DG.DataFrame.fromProperties(f.inputs.concat(f.outputs), 5);
    let grid = DG.Grid.create(table);
    let calculating: boolean = false;

    function refreshRow(row: number) {
      for (let i = 0; i < f.inputs.length; i++)
        if (table.columns.byIndex(i).isNone(row))
          return;

      let params: any = {};
      for (let p of f.inputs)
        params[p.name] = table.col(p.name)!.get(row);

      calculating = true;
      let call = f.prepare(params);
      call.call().then(c => {
        for (let p of f.outputs)
          table.col(p.name)!.set(row, call.outputs[p.name]);
        calculating = false;
      });
    }

    grid.onCellValueEdited.subscribe(gc => refreshRow(gc.tableRow?.idx!));
    gridHost.appendChild(grid.root);
  }

  init();
  functionName.onChanged(() => {
    DG.Func
      .findAll({name: functionName.value})
      .then(funcs => {
        if (funcs.length == 1) {
          f = funcs[0];
          init();
        }
      });
    // ;
    // if (f !== null)
    //   init();
  });
  return view;
}