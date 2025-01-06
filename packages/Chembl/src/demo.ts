import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function _demoDatabasesChembl(): Promise<void> {
  const query = `--name: compound activity details for target 
--connection: Chembl
--input: string target_name = "Acetylcholinesterase" {choices: Query("SELECT distinct pref_name from target_dictionary limit 300 offset 309;")}
--input: string target_id = '93' {choices: Query("SELECT distinct tid::int from target_dictionary where pref_name = @target_name;")}
--input: string substructure = "NC1=CC(=O)c2ccccc2C1=O" {semType: Substructure}
--input: string activity_type = "IC50"
  
SELECT canonical_smiles, description, standard_inchi, t.target_type, c.molregno, a.chembl_id, a.assay_id FROM assays a  
JOIN target_dictionary t on a.tid = t.tid 
JOIN activities act on a.assay_id = act.assay_id
JOIN compound_structures c on act.molregno = c.molregno
WHERE t.tid = CAST(@target_id as integer)
AND c.canonical_smiles::mol @>@substructure::qmol
AND act.type = @activity_type
LIMIT 50
--end`;

  const connection: DG.DataConnection = await grok.functions.eval('Chembl:Chembl');

  const dBQuery = connection!.query('', query);
  const funccall = dBQuery.prepare();
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
