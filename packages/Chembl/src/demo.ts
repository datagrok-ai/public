import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function _demoDatabasesChembl(): Promise<void> {
  const query: DG.DataQuery = await grok.functions.eval('Chembl:FracClassificationWithSubstructure');
  const funccall: DG.FuncCall = query.prepare();
  const editor = await funccall.getEditor();

  const queryPanel = ui.input.textArea('', {value: query.query});
  queryPanel.input.style.width = '100%';
  queryPanel.input.style.minHeight = '350px';
  queryPanel.input.setAttribute('readonly', 'true');
  const gridDiv = ui.div('', {style: {position: 'relative', height: '100%'}});

  const runQuery = async () => {
    ui.setUpdateIndicator(gridDiv, true);
    await funccall.call();
    const data: DG.DataFrame = funccall.getOutputParamValue();
    if (data) {
      await grok.data.detectSemanticTypes(data);
      const grid = data.plot.grid().root;
      grid.style.width = '100%';
      grid.style.height = '100%';
      ui.empty(gridDiv);
      gridDiv.append(grid);
      ui.setUpdateIndicator(gridDiv, false);
    }
  };

  const tabControl = ui.tabControl({
    'Query Input Form': ui.divV([editor]),
    'Query SQL': ui.div(queryPanel.root),
  });
  tabControl.root.style.width = '100%';
  tabControl.root.style.height = '375px';

  const totalDiv = ui.divV([
    tabControl.root,
    gridDiv,
  ], {style: {height: '100%', width: '100%'}});

  const view = DG.View.create();
  view.name = 'Chemical Databases';
  view.root.append(totalDiv);

  for (const p of funccall.inputParams.values())
    view.subs.push(p.onChanged.subscribe(() => runQuery()));

  grok.shell.addView(view);
  await runQuery();
}
