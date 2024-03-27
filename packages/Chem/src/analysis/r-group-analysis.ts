import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {findRGroups, findRGroupsWithCore} from '../scripts-api';
import {convertMolNotation, getRdKitModule} from '../package';
import {getMCS} from '../utils/most-common-subs';
import {IRGroupAnalysisResult} from '../rdkit-service/rdkit-service-worker-substructure';
import {getRdKitService} from '../utils/chem-common-rdkit';
import { MolfileHandler } from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { _convertMolNotation } from '../utils/convert-notation-utils';
import { SCAFFOLD_COL, SUBSTRUCT_COL } from '../constants';
import { getUncommonAtomsAndBonds } from '../utils/chem-common-rdkit';
import { delay } from '@datagrok-libraries/utils/src/test';

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

let latestAnalysisCols: {[key: string]: string []} = {};
let latestTrellisPlot: {[key: string]: DG.Viewer | null} = {};

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
  const replaceLatest = ui.boolInput('Replace latest', true);
  const undoLatestAnalysis = ui.icons.undo(() => {
    removeLatestAnalysis(col);
  }, 'Undo latest analysis');
  undoLatestAnalysis.classList.add('chem-rgroup-undo-icon');

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
      ui.divH([visualAnalysisCheck.root, replaceLatest.root, undoLatestAnalysis])
    ]))
    .onOK(async () => {
      if (replaceLatest.value) {
        removeLatestAnalysis(col);
        await delay(50);
      }
      const getPrefixIdx = (colPrefix: string) => {
        let prefixIdx = 0;
        col = col.dataFrame.columns.byName(columnInput.value!);
        const re = new RegExp(`^${colPrefix}$`, 'i');
        if (col.dataFrame.columns.names().filter(((it) => it.match(re))).length) {
          prefixIdx++;
          const maxPrefixIdx = 100;
          for (let i = 0; i < maxPrefixIdx; i++) {
            const reIdx = new RegExp(`^${colPrefix}_${prefixIdx}$`, 'i');
            if (!col.dataFrame.columns.names().filter(((it) => it.match(reIdx))).length)
              break;
            prefixIdx++;
          }
          if (prefixIdx - 1 === maxPrefixIdx) {
            grok.shell.error('Table contains columns named \'R[number]\', please change column prefix');
            return null;
          }
        }
        return prefixIdx;
      }
      const rGroupPrefixRe = `${columnPrefixInput.value}\\d+`;
      const corePrefixRe = `Core`;
      const rGroupPrefixIdx = getPrefixIdx(rGroupPrefixRe);
      const corePrefixIdx = getPrefixIdx(corePrefixRe);
      if (rGroupPrefixIdx === null || corePrefixIdx === null)
        return;
      let core = sketcher.getMolFile();
      if (!core) {
        grok.shell.error('No core was provided');
        return;
      }
      let progressBar;
      try {
        const onlyMatchAtRGroups =
          !!MolfileHandler.getInstance(core).atomTypes.filter((it) => it.startsWith('R')).length &&
          core.includes('M  RGP');
        if (!onlyMatchAtRGroups)
          core = convertMolNotation(core, grok.chem.Notation.MolBlock, grok.chem.Notation.Smarts);
        progressBar = DG.TaskBarProgressIndicator.create(`RGroup analysis running...`);
        //const res = await rGroupsPython(col, core, columnPrefixInput.value, true, onlyMatchAtRGroups);

        const rGroupOptions = {
          chunkSize: rGroupChunkSize.value!.toString(),
          matchingStrategy: rGroupMatchingStrategy.value!,
          alignment: rGroupAlignment.value!,
        };
        const res = await rGroupsMinilib(col, core, rGroupOptions);

        const module = getRdKitModule();
        if (res.length) {
          //unmatched are those items for which all R group cols are empty
          const unmatchedItems = new Uint8Array(res[0].length).fill(0);
          latestAnalysisCols[col.dataFrame.name] = [];
          for (const resCol of res) {
            const molsArray = new Array<string>(resCol.length);
            for (let i = 0; i < resCol.length; i++) {
              const molStr = resCol.get(i);
              if (resCol.name !== 'Core' && !molStr) {
                unmatchedItems[i] += 1;
              } else {
                let mol: RDMol | null = null;
                try {
                  mol = module.get_mol(molStr!);
                  if (mol)
                    molsArray[i] = mol.get_molblock().replace('ISO', 'RGP');
                } catch (e) {
                  //do nothing here, molsArray[i] is empty for invalid molecules
                } finally {
                  mol?.delete();
                }
              }
            }
            let rColName = '';
            if (resCol.name === 'Core') {
              rColName = corePrefixIdx ? `${resCol.name}_${corePrefixIdx}` : resCol.name;
              col.temp[SCAFFOLD_COL] = rColName;
            }
            else
              rColName = rGroupPrefixIdx ? `${resCol.name}_${rGroupPrefixIdx}` : resCol.name;
            resCol.name = rColName;
            const rCol = DG.Column.fromStrings(rColName, molsArray);
            rCol.semType = DG.SEMTYPE.MOLECULE;
            rCol.setTag(DG.TAGS.UNITS, DG.chem.Notation.MolBlock);
            col.dataFrame.columns.add(rCol);
            latestAnalysisCols[col.dataFrame.name].push(rColName);
          }
          //create column for highlight of uncommon structure
          const highlightColName = `r-groups-highlight_${rGroupPrefixIdx}`;
          let coreMol: RDMol | null = null;
          let coreMolWithoutRGroups: RDMol | null = null;
          try {
            //const coreWithoutRGroups = core.replaceAll('R#', '*');
            coreMol = module.get_mol(core);
            const smiles = coreMol.get_smiles();
            //remove r groups from smiles string to make highlight more precise
            const newSmiles = smiles.replaceAll(/(\[\d+\*\])/g, '').replaceAll('()', '');
            coreMolWithoutRGroups = module.get_mol(newSmiles);
            const substructCol = DG.Column.fromType('object', highlightColName, col.dataFrame.rowCount)
              .init((i) => getUncommonAtomsAndBonds(col.get(i), coreMolWithoutRGroups, module, '#bc131f', true));
            col.dataFrame.columns.add(substructCol);
            latestAnalysisCols[col.dataFrame.name].push(highlightColName);
            col.temp[SUBSTRUCT_COL] = highlightColName;
          } finally {
            coreMol?.delete();
            coreMolWithoutRGroups?.delete();
          }
          //create boolean column for match/non match
          const isHitCol = DG.Column.bool(`${rGroupPrefixIdx ? `isHit_${rGroupPrefixIdx}` : `isHit`}`, res[0].length)
            .init((i) => unmatchedItems[i] !== res.length - 1);
          col.dataFrame.columns.add(isHitCol);
          latestAnalysisCols[col.dataFrame.name].push(isHitCol.name);
          //filter out unmatched values
          const filterUnmatched = DG.BitSet.create(res[0].length).init((i) => isHitCol.get(i));
          col.dataFrame.filter.copyFrom(filterUnmatched);
          const view = grok.shell.getTableView(col.dataFrame.name);
          //make highlight column invisible
          view.grid.col(highlightColName)!.visible = false;
          if (visualAnalysisCheck.value && view) {
            if (res.length < 3)
              grok.shell.warning(`Not enough R group columns to create trellis plot`);
            else
              latestTrellisPlot[col.dataFrame.name] = view.trellisPlot({
                xColumnNames: [res[1].name], // column 0 is Core column
                yColumnNames: [res[2].name],
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
  const mod = getRdKitModule();
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

export async function rGroupsPython(col: DG.Column<string>, core: string, prefix: string, withCore: boolean,
  onlyMatchAtRGroups: boolean):
  Promise<DG.Column<string>[]> {
  const resCols = [];
  const res = withCore ? await findRGroupsWithCore(col.name, col.dataFrame, core, onlyMatchAtRGroups):
    await findRGroups(col.name, col.dataFrame, core, prefix);
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
function removeLatestAnalysis(col: DG.Column) {
  if(latestTrellisPlot[col.dataFrame.name] && latestTrellisPlot[col.dataFrame.name]!.dataFrame) 
    latestTrellisPlot[col.dataFrame.name]!.close();
  delete latestTrellisPlot[col.dataFrame.name];
  if(latestAnalysisCols[col.dataFrame.name])
    latestAnalysisCols[col.dataFrame.name]!.forEach((colName: string) => col.dataFrame.columns.remove(colName));
  delete latestAnalysisCols[col.dataFrame.name];
}
