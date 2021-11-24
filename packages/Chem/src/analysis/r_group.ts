import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { mcsgetter } from '../scripts-api';
import { getRGroups } from '../chem_rgroup_analysis';

export function rGroupAnalysis(col: DG.Column) {
  let sketcherSmile = '';

  let sketcher = new grok.chem.Sketcher();
  let columnPrefixInput = ui.stringInput('Column prefix', 'R');
  let visualAnalysisCheck = ui.boolInput('Visual analysis', true);

  let mcsButton = ui.button('MCS', async () => {
    let smiles: string = await mcsgetter(col.name, col.dataFrame);
    sketcher.setSmiles(smiles);
    sketcherSmile = smiles;
  });
  ui.tooltip.bind(mcsButton, "Most Common Substructure");
  let mcsButtonHost = ui.div([mcsButton]);
  mcsButtonHost.style.display = 'flex';
  mcsButtonHost.style.justifyContent = 'center';

  let dlg = ui.dialog({
    title: 'R-Groups Analysis',
    helpUrl: '/help/domains/chem/cheminformatics.md#r-group-analysis'
    })
    .add(ui.div([
      sketcher,
      mcsButtonHost,
      columnPrefixInput,
      visualAnalysisCheck
    ]))
    .onOK(async () => {
      let res = await getRGroups(col, sketcherSmile, columnPrefixInput.value);
      for (let resCol of res.columns) {
        resCol.semType = DG.SEMTYPE.MOLECULE;
        col.dataFrame.columns.add(resCol);
      }
      if (res.columns.length == 0)
        grok.shell.error("None R-Groups were found");
      let view = grok.shell.getTableView(col.dataFrame.name);
      if (visualAnalysisCheck.value && view) {
        view.trellisPlot({
          xColumnNames: [res.columns[0].name],
          yColumnNames: [res.columns[1].name]});
      }
    });
  dlg.show();
  dlg.initDefaultHistory();
}