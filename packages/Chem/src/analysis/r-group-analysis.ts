import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {findRGroups} from '../scripts-api';
import {getRdKitModule} from '../package';
import {getMCS} from '../utils/most-common-subs';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';


export function convertToRDKit(smiles: string): string {
  const regexConv: RegExp = /(\[)(R)(\d+)(\])/g;
  const match = regexConv.exec(smiles);
  if (match !== null)
    smiles = smiles.replace(regexConv, `${match[1]}*:${match[3]}${match[4]}`);

  return smiles;
}

/**
 * Opens a dialog with the sketcher that allows to sketch a core (or perform MCS),
 * and initiate the R-Group Analysis for the specified column with molecules.
 * @param {DG.Column} col Column for which to perform R-Group Analysis with specified core
 */
export function rGroupAnalysis(col: DG.Column): void {
  const sketcher = new DG.chem.Sketcher();
  const columnPrefixInput = ui.stringInput('Column prefix', 'R');
  const visualAnalysisCheck = ui.boolInput('Visual analysis', true);
  const exactAtomsCheck = ui.boolInput('MCS exact atoms', true);
  const exactBondsCheck = ui.boolInput('MCS exact bonds', true);

  const molColNames = col.dataFrame.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).map((c) => c.name);
  const columnInput = ui.choiceInput('Molecules', col.name, molColNames);

  const mcsButton = ui.button('MCS', async () => {
    ui.setUpdateIndicator(sketcher.root, true);
    try {
      const molCol = col.dataFrame.columns.byName(columnInput.value!);
      //TODO: implements mcs using web worker
      const mcsSmarts = await getMCS(molCol, exactAtomsCheck.value!, exactBondsCheck.value!);
      if (mcsSmarts !== null) {
        ui.setUpdateIndicator(sketcher.root, false);
        sketcher.setSmarts(mcsSmarts);
      }
    } catch (e: any) {
      grok.shell.error(e);
      dlg.close();
    }
  });
  ui.tooltip.bind(mcsButton, 'Calculate Most Common Substructure');
  const mcsButtonHost = ui.div([mcsButton], 'chem-mcs-button-host');
  mcsButton.classList.add('chem-mcs-button');

  const dlg = ui.dialog({
    title: 'R-Groups Analysis',
    helpUrl: '/help/datagrok/solutions/domains/chem/#r-groups-analysis',
  })
    .add(ui.div([
      sketcher,
      ui.divH([
        ui.divV([
          exactAtomsCheck.root, exactBondsCheck.root,
        ]),
        mcsButtonHost,
      ]),
      columnInput,
      columnPrefixInput,
      visualAnalysisCheck.root,
    ]))
    .onOK(async () => {
      col = col.dataFrame.columns.byName(columnInput.value!);
      const re = new RegExp(`^${columnPrefixInput.value}\\d+$`, 'i');
      if (col.dataFrame.columns.names().filter(((it) => it.match(re))).length) {
        grok.shell.error('Table contains columns named \'R[number]\', please change column prefix');
        return;
      }
      const core = await sketcher.getSmarts();
      if (!core) {
        grok.shell.error('No core was provided');
        return;
      }
      let progressBar;
      try {
        progressBar = DG.TaskBarProgressIndicator.create(`RGroup analysis running...`);
        const res = await findRGroups(col.name, col.dataFrame, core, columnPrefixInput.value);
        const module = getRdKitModule();
        if (res.rowCount) {
          for (const resCol of res.columns) {
            const molsArray = new Array<string>(resCol.length);
            for (let i = 0; i < resCol.length; i++) {
              const molStr = resCol.get(i);
              let mol: RDMol | null = null;
              try {
                mol = module.get_mol(molStr);
                molsArray[i] = mol.get_molblock().replace('ISO', 'RGP');
              } catch (e) {
                //do nothing here, molsArray[i] is empty for invalid molecules
              } finally {
                mol?.delete();
              }
            }
            const rCol = DG.Column.fromStrings(resCol.name, molsArray);

            rCol.semType = DG.SEMTYPE.MOLECULE;
            rCol.setTag(DG.TAGS.UNITS, DG.chem.Notation.MolBlock);
            col.dataFrame.columns.add(rCol);
          }

          const view = grok.shell.getTableView(col.dataFrame.name);
          if (visualAnalysisCheck.value && view) {
            view.trellisPlot({
              xColumnNames: [res.columns.byIndex(0).name],
              yColumnNames: [res.columns.byIndex(1).name],
            });
          }
        } else
          grok.shell.error('None R-Groups were found');
        progressBar.close();
      } catch (e: any) {
        grok.shell.error(e);
        dlg.close();
        progressBar?.close();
      }
    });
  dlg.show();
  dlg.initDefaultHistory();
}
