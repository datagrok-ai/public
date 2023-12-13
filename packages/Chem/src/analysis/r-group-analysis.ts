import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {findRGroups} from '../scripts-api';
import {getRdKitModule} from '../package';
import {getMCS} from '../utils/most-common-subs';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IRGroupAnalysisResult} from '../rdkit-service/rdkit-service-worker-substructure';
import {getRdKitService} from '../utils/chem-common-rdkit';

const enum RGroupMatchingStrategy {
  Greedy = 'Greedy',
  GreedyChunks = 'GreedyChunks',
  Exhaustive = 'Exhaustive',
  NoSymmetrization = 'NoSymmetrization',
  GA = 'GA'
};

const enum RGroupAlignment {
  None = 'None',
  NoAlignment = 'NoAlignment',
  MCS = 'MCS'
};

const matchingStrategies: RGroupMatchingStrategy[] = [
  RGroupMatchingStrategy.Greedy,
  RGroupMatchingStrategy.GreedyChunks,
  RGroupMatchingStrategy.Exhaustive,
  RGroupMatchingStrategy.GA,
  RGroupMatchingStrategy.NoSymmetrization,
];

const alignments: RGroupAlignment[] = [RGroupAlignment.MCS, RGroupAlignment.NoAlignment, RGroupAlignment.None];


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

  const rGroupMatchingStrategy = ui.choiceInput('Matching strategy', matchingStrategies[0], matchingStrategies);
  const rGroupAlignment = ui.choiceInput('Alignment', alignments[0], alignments);
  const rGroupChunkSize = ui.intInput('Chunk size', 5);

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


  const rGroupsSettings = ui.icons.settings(() => {
    rGroupSettinsOpened = !rGroupSettinsOpened;
    if (!rGroupSettinsOpened)
      ui.empty(rGroupSettingsDiv);
    else {
      rGroupSettingsDiv.append(rGroupMatchingStrategy.root);
      rGroupSettingsDiv.append(rGroupAlignment.root);
      rGroupSettingsDiv.append(rGroupChunkSize.root);
    }
  }, 'R group analysis settins');
  rGroupsSettings.classList.add('r-group-settings-icon');
  const rGroupAnalysisDiv = ui.divH([
    ui.divText('R Group Analysis', 'r-group-settins-div'),
    rGroupsSettings,
  ]);
  const rGroupSettingsDiv = ui.div([], 'rgroup-mcs-settings');
  let rGroupSettinsOpened = false;


  const dlg = ui.dialog({
    title: 'R-Groups Analysis',
    helpUrl: '/help/datagrok/solutions/domains/chem/chem.md#r-groups-analysis',
  })
    .add(ui.div([
      sketcher,
      ui.divV([
        ui.divH([mcsButton, exactAtomsCheck.root, exactBondsCheck.root]),
        rGroupAnalysisDiv,
        rGroupSettingsDiv,
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
        //const result = await rGroupsPython(col, core, columnPrefixInput.value);
        const rGroupOptions = {
          chunkSize: rGroupChunkSize.value!.toString(),
          matchingStrategy: rGroupMatchingStrategy.value!,
          alignment: rGroupAlignment.value!,
        };
        const result = await rGroupsMinilib(col, core, rGroupOptions);

        if (result.length) {
          result.forEach((rCol) => col.dataFrame.columns.add(rCol));


          const view = grok.shell.getTableView(col.dataFrame.name);
          if (visualAnalysisCheck.value && view) {
            view.trellisPlot({
              xColumnNames: [result[0].name],
              yColumnNames: [result[1].name],
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

export async function rGroupsMinilib(molecules: DG.Column<string>, coreMolecule: string, options?:
    {[key: string]: string}): Promise<DG.Column<string>[]> {
  const res: IRGroupAnalysisResult =
    await (await getRdKitService())
      .getRGroups(molecules.toList(), coreMolecule, options ? JSON.stringify(options) : '');
  const resCols: DG.Column<string>[] = [];
  for (let i = 0; i < res.colNames.length; i++) {
    const col = DG.Column.string(res.colNames[i], molecules.length).init((j) => res.smiles[i][j]);
    col.semType = DG.SEMTYPE.MOLECULE;
    resCols.push(col);
  }
  return resCols;
}

export async function rGroupsPython(col: DG.Column<string>, core: string, prefix: string):
  Promise<DG.Column<string>[]> {
  const resCols = [];
  const res = await findRGroups(col.name, col.dataFrame, core, prefix);
  const module = getRdKitModule();
  if (res.rowCount) {
    for (const resCol of res.columns) {
      const molsArray = new Array<string>(resCol.length);
      for (let i = 0; i < resCol.length; i++) {
        const molStr = resCol.get(i);
        try {
          const mol = module.get_mol(molStr);
          molsArray[i] = mol.get_molblock().replace('ISO', 'RGP');
          mol.delete();
        } catch (e) {
          console.warn(`RGroupAnalysisWarning: skipping invalid molecule '${molStr}' at index ${i}`);
        }
      }
      const rCol = DG.Column.fromStrings(resCol.name, molsArray);
      rCol.semType = DG.SEMTYPE.MOLECULE;
      rCol.setTag(DG.TAGS.UNITS, DG.chem.Notation.MolBlock);
      resCols.push(rCol);
    }
  }
  return resCols;
}
