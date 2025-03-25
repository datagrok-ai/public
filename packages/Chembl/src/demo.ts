import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function _demoDatabasesChembl(): Promise<void> {
  const query = await grok.functions.eval('Chembl:FracClassificationWithSubstructure');
  const funccall = query.prepare();
  const editor = await funccall.getEditor();
  const runButton = ui.bigButton('RUN', async () => {
    await runQuery();
  });

  const runQuery = async () => {
    ui.setUpdateIndicator(gridDiv, true);
    await funccall.call();
    const data: DG.DataFrame = funccall.getOutputParamValue();
    await grok.data.detectSemanticTypes(data);
    const grid = data.plot.grid().root;
    grid.style.width = '100%';
    grid.style.height = '100%';
    ui.empty(gridDiv);
    gridDiv.append(grid);
    ui.setUpdateIndicator(gridDiv, false);
  };

  runButton.style.width = '150px';
  runButton.style.marginLeft = '80px';

  const queryPanel = ui.input.textArea('', {value: query});
  queryPanel.input.style.width = '100%';
  queryPanel.input.style.minHeight = '350px';
  const gridDiv = ui.div('', {style: {position: 'relative', height: '100%'}});

  const tabControl = ui.tabControl({
    'Query Input Form': ui.divV([
      editor,
      runButton,
    ]),
    'Query SQL': queryPanel,
  });
  tabControl.root.style.width = '100%';
  tabControl.root.style.height = '310px';

  const totalDiv = ui.divV([
    tabControl.root,
    gridDiv,
  ], {style: {height: '100%', width: '100%'}});

  const view = grok.shell.addView(DG.View.create());
  view.name = 'Chemical Databases';
  view.root.append(totalDiv);
  setTimeout(() => runQuery(), 0);
}
