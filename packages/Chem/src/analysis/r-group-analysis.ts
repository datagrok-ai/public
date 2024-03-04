import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {findRGroups, findRGroupsWithCore} from '../scripts-api';
import {convertMolNotation, getRdKitModule} from '../package';
import {getMCS} from '../utils/most-common-subs';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { _convertMolNotation } from '../utils/convert-notation-utils';
import { SCAFFOLD_COL, SUBSTRUCT_COL } from '../constants';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import { getUncommonAtomsAndBonds } from '../utils/chem-common-rdkit';


let latestAnalysisCols: string[] = [];
let latestTrellisPlot: DG.Viewer | null = null;

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
    helpUrl: '/help/datagrok/solutions/domains/chem/chem.md#r-groups-analysis',
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
      ui.divH([visualAnalysisCheck.root, replaceLatest.root, undoLatestAnalysis])
    ]))
    .onOK(async () => {
      if (replaceLatest.value)
        removeLatestAnalysis(col);
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
        const res = await findRGroupsWithCore(col.name, col.dataFrame, core, onlyMatchAtRGroups);
        const module = getRdKitModule();
        if (res.rowCount) {
          //unmatched are those items for which all R group cols are empty
          const unmatchedItems = new Uint8Array(res.rowCount).fill(0);
          for (const resCol of res.columns) {
            const molsArray = new Array<string>(resCol.length);
            for (let i = 0; i < resCol.length; i++) {
              const molStr = resCol.get(i);
              if (resCol.name !== 'Core' && !molStr) {
                unmatchedItems[i] += 1;
              } else {
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
            latestAnalysisCols.push(rColName);
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
            latestAnalysisCols.push(highlightColName);
            col.temp[SUBSTRUCT_COL] = highlightColName;
          } finally {
            coreMol?.delete();
            coreMolWithoutRGroups?.delete();
          }
          //create boolean column for match/non match
          const isHitCol = DG.Column.bool(`${rGroupPrefixIdx ? `isHit_${rGroupPrefixIdx}` : `isHit`}`, res.rowCount)
            .init((i) => unmatchedItems[i] !== res.columns.length - 1);
          col.dataFrame.columns.add(isHitCol);
          latestAnalysisCols.push(isHitCol.name);
          //filter out unmatched values
          const filterUnmatched = DG.BitSet.create(res.rowCount).init((i) => isHitCol.get(i));
          col.dataFrame.filter.copyFrom(filterUnmatched);
          const view = grok.shell.getTableView(col.dataFrame.name);
          //make highlight column invisible
          view.grid.col(highlightColName)!.visible = false;
          if (visualAnalysisCheck.value && view) {
            latestTrellisPlot = view.trellisPlot({
              xColumnNames: [res.columns.byIndex(1).name], // column 0 is Core column
              yColumnNames: [res.columns.byIndex(2).name],
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

function removeLatestAnalysis(col: DG.Column) {
  latestTrellisPlot?.close();
  latestTrellisPlot = null;
  latestAnalysisCols.forEach((colName: string) => col.dataFrame.columns.remove(colName));
  latestAnalysisCols = [];
}