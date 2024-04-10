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
import { delay } from '@datagrok-libraries/utils/src/test';
import { hexToPercentRgb } from '../utils/chem-common';
import { ISubstruct } from '../rendering/rdkit-cell-renderer';

export type RGroupParams = {
  molColName: string,
  core: string,
  rGroupName: string,
  rGroupChunkSize: string,
  rGroupMatchingStrategy: string,
  rGroupAlignment: string,
  visualAnalysis: boolean,
}

export type RGroupDecompRes = {
  xAxisColName: string,
  yAxisColName: string,
  highlightColName?: string,
}

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

type RGroupsRes = {
  rGroups: DG.Column<string>[];
  highlightCol?: DG.Column<object>;
}

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

  //General fields 
  const molColNames = col.dataFrame.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).map((c) => c.name);
  const columnInput = ui.choiceInput('Molecules', col.name, molColNames);
  const columnPrefixInput = ui.stringInput('Column prefix', 'R');
  const visualAnalysisCheck = ui.boolInput('Visual analysis', true);
  const replaceLatest = ui.boolInput('Replace latest', true);
  replaceLatest.root.classList.add('chem-rgroup-replace-latest');

  //MCS fields
  const mcsButton = ui.button('MCS', async () => {
    ui.setUpdateIndicator(sketcher.root, true);
    try {
      const molCol = col.dataFrame.columns.byName(columnInput.value!);
      //TODO: implements mcs using web worker
      const mcsSmarts = await getMCS(molCol, mcsExactAtomsCheck.value!, mcsExactBondsCheck.value!);
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
  const mcsExactAtomsCheck = ui.boolInput('Exact atoms', true);
  const mcsExactBondsCheck = ui.boolInput('Exact bonds', true);
  mcsExactBondsCheck.captionLabel.classList.add('chem-mcs-settings-label');

  //R groups fields
  const rGroupMatchingStrategy = ui.choiceInput('Matching strategy', matchingStrategies[0], matchingStrategies);
  const rGroupAlignment = ui.choiceInput('Alignment', alignments[0], alignments);
  const rGroupChunkSize = ui.intInput('Chunk size', 5);

  //settings button to adjust mcs and r-groups settings
  const mcsAndrGroupsSettingsIcon = ui.iconFA('cog', () => {
    rGroupSettinsOpened = !rGroupSettinsOpened;
    if (!rGroupSettinsOpened)
      ui.empty(mcsAndRGroupSettingsDiv);
    else {
      mcsAndRGroupSettingsDiv.append(rGroupMatchingStrategy.root);
      mcsAndRGroupSettingsDiv.append(rGroupAlignment.root);
      mcsAndRGroupSettingsDiv.append(rGroupChunkSize.root);
    }
  }, 'R group analysis settings');
  mcsAndrGroupsSettingsIcon.classList.add('chem-rgroup-settings-icon');
  const mcsAndRGroupSettingsDiv = ui.inputs([], 'chem-rgroup-settings-div');
  let rGroupSettinsOpened = false;


  const dlg = ui.dialog({
    title: 'R-Groups Analysis',
    helpUrl: '/help/datagrok/solutions/domains/chem/chem.md#r-groups-analysis',
  })
    .add(ui.divV([
      sketcher,
      ui.divH([mcsButton, mcsExactAtomsCheck.root, mcsExactBondsCheck.root], {style: {paddingLeft: '108px'}}),
      columnInput,
      columnPrefixInput,
      visualAnalysisCheck.root,
      mcsAndrGroupsSettingsIcon,
      latestAnalysisCols[col.dataFrame.name]?.length ? replaceLatest.root : null,
      mcsAndRGroupSettingsDiv
    ]))
    .onOK(async () => {
      try {
        if (replaceLatest.value) {
          removeLatestAnalysis(col);
          await delay(50);
        }
        const smarts = await sketcher.getSmarts();
        const funcCall = await DG.Func.find({ name: 'rGroupDecomposition' })[0].prepare({
          df: col.dataFrame,
          molColName: columnInput.value!,
          core: smarts,
          rGroupName: columnPrefixInput.value,
          rGroupChunkSize: rGroupChunkSize.value!.toString(),
          rGroupMatchingStrategy: rGroupMatchingStrategy.value!,
          rGroupAlignment: rGroupAlignment.value!,
          visualAnalysis: visualAnalysisCheck.value!,
        }).call(undefined, undefined, { processed: false });
        const res: RGroupDecompRes = funcCall.getOutputParamValue();
        if (res) {
          const view = grok.shell.getTableView(col.dataFrame.name);
          //make highlight column invisible
          if (res.highlightColName)
            view.grid.col(res.highlightColName)!.visible = false;
          if (visualAnalysisCheck.value! && view) {
            if (!res.yAxisColName || !res.xAxisColName)
              grok.shell.warning(`Not enough R group columns to create trellis plot`);
            else
              latestTrellisPlot[col.dataFrame.name] = view.trellisPlot({
                xColumnNames: [res.xAxisColName],
                yColumnNames: [res.yAxisColName],
              });
          }
        }
      } catch (e: any) {
        grok.shell.error(e.message);
      }
    });
  dlg.show();
  dlg.initDefaultHistory();
}

export async function rGroupDecomp(col: DG.Column, params: RGroupParams): Promise<RGroupDecompRes | undefined> {
  const getPrefixIdx = (colPrefix: string) => {
    let prefixIdx = 0;
    col = col.dataFrame.columns.byName(params.molColName);
    const re = new RegExp(`^${colPrefix}$`, 'i');
    if (col.dataFrame.columns.names().filter(((it) => it.match(re))).length) {
      prefixIdx++;
      const maxPrefixIdx = 2;
      for (let i = 0; i < maxPrefixIdx; i++) {
        const reIdx = new RegExp(`^${colPrefix}_${prefixIdx}$`, 'i');
        if (!col.dataFrame.columns.names().filter(((it) => it.match(reIdx))).length)
          break;
        prefixIdx++;
      }
      if (prefixIdx - 1 === maxPrefixIdx) {
        throw new Error('Table contains columns named \'R[number]\', please change column prefix');
      }
    }
    return prefixIdx;
  }
  const rGroupPrefixRe = `${params.rGroupName}\\d+`;
  const corePrefixRe = `Core`;
  const rGroupPrefixIdx = getPrefixIdx(rGroupPrefixRe);
  const corePrefixIdx = getPrefixIdx(corePrefixRe);
  let core = params.core;
  if (!core) {
    throw new Error('No core was provided');
  }
  let progressBar;
  try {
    //using onlyMatchAtRGroups param in case user has manually defined enumerated R groups
    core = convertMolNotation(core, grok.chem.Notation.Smarts, grok.chem.Notation.MolBlock);
    const onlyMatchAtRGroups =
      !!MolfileHandler.getInstance(core).atomTypes.filter((it) => it.startsWith('R')).length &&
      core.includes('M  RGP');
    const coreIsQMol = core.includes('M  ALS') || core.includes('M  RAD');
    if (coreIsQMol)
      core = convertMolNotation(core, grok.chem.Notation.MolBlock, grok.chem.Notation.Smarts);
    progressBar = DG.TaskBarProgressIndicator.create(`RGroup analysis running...`);
    //const res = await rGroupsPython(col, core, columnPrefixInput.value, true, onlyMatchAtRGroups);

    const rGroupOptions = {
      chunkSize: params.rGroupChunkSize,
      matchingStrategy: params.rGroupMatchingStrategy,
      alignment: params.rGroupAlignment,
      includeTargetMolInResults: true,
      onlyMatchAtRGroups: onlyMatchAtRGroups
    };
    const { rGroups, highlightCol } = await rGroupsMinilib(col, core, coreIsQMol, rGroupPrefixIdx, rGroupOptions);
    const module = getRdKitModule();
    if (rGroups.length) {
      //unmatched are those items for which all R group cols are empty
      const unmatchedItems = new Uint8Array(rGroups[0].length).fill(0);
      latestAnalysisCols[col.dataFrame.name] = [];
      for (const resCol of rGroups) {
        const molsArray = new Array<string>(resCol.length);
        for (let i = 0; i < resCol.length; i++) {
          const molStr = resCol.get(i);
          if (resCol.name !== 'Core' && !molStr) {
            unmatchedItems[i] += 1;
          } else {
            let mol: RDMol | null = null;
            try {
              mol = resCol.name === 'Core' && coreIsQMol ? module.get_qmol(molStr!) : module.get_mol(molStr!);
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
      //create column for r groups highlight
      if (highlightCol) {
        col.dataFrame.columns.add(highlightCol);
        latestAnalysisCols[col.dataFrame.name].push(highlightCol.name);
        col.temp[SUBSTRUCT_COL] = highlightCol.name;
      }
      //create boolean column for match/non match
      const isHitCol = DG.Column.bool(`${rGroupPrefixIdx ? `isHit_${rGroupPrefixIdx}` : `isHit`}`, rGroups[0].length)
        .init((i) => unmatchedItems[i] !== rGroups.length - 1);
      col.dataFrame.columns.add(isHitCol);
      latestAnalysisCols[col.dataFrame.name].push(isHitCol.name);
      //filter out unmatched values
      const filterUnmatched = DG.BitSet.create(rGroups[0].length).init((i) => isHitCol.get(i));
      col.dataFrame.filter.copyFrom(filterUnmatched);
    } else
      grok.shell.warning('None R-Groups were found');
    progressBar.close();
    return {
      xAxisColName: rGroups.length > 1 ? rGroups[1].name : '',  //rGroups[0] column is Core column
      yAxisColName: rGroups.length > 2 ? rGroups[2].name : '',
      highlightColName: highlightCol?.name
    };
  } catch (e: any) {
    grok.shell.error(e);
    progressBar?.close();
  }
}

export async function rGroupsMinilib(molecules: DG.Column<string>, coreMolecule: string,
  coreIsQMol: boolean, rGroupPrefixIdx: number, options?:
    { [key: string]: string | boolean }): Promise<RGroupsRes> {
  if (!coreMolecule)
      throw new Error('No core was provided');
  const res: IRGroupAnalysisResult =
    await (await getRdKitService())
      .getRGroups(molecules.toList(), coreMolecule, coreIsQMol, options ? JSON.stringify(options) : '');
  const resCols: DG.Column<string>[] = [];
  for (let i = 0; i < res.colNames.length; i++) {
    const col = DG.Column.string(res.colNames[i], molecules.length).init((j) => res.smiles[i][j]);
    col.semType = DG.SEMTYPE.MOLECULE;
    resCols.push(col);
  }
  //creating highlight column
  const highlightColName = `r-groups-highlight_${rGroupPrefixIdx}`;
  const rGroupsNum = res.atomsToHighLight.length;
  const colors = Array<Array<number>>(rGroupsNum);
  for (let r = 0; r < rGroupsNum; r++)
    colors[r] = hexToPercentRgb(DG.Color.toHtml(DG.Color.getCategoricalColor(r)))!;

  const substructCol = DG.Column.fromType('object', highlightColName, molecules.dataFrame.rowCount)
    .init((i) => {
      const substr: ISubstruct = {
        atoms: [],
        bonds: [],
        highlightAtomColors: {},
        highlightBondColors: {},
      };
      for (let j = 0; j < rGroupsNum; j++) {
        if(res.atomsToHighLight[j][i]) {
          res.atomsToHighLight[j][i].forEach((atom) => substr.highlightAtomColors![atom] = colors[j]);
          substr.atoms = substr.atoms!.concat(Array.from(res.atomsToHighLight[j][i]));
        }
        if(res.bondsToHighLight[j][i]) {
          res.bondsToHighLight[j][i].forEach((bond) => substr.highlightBondColors![bond] = colors[j]);
          substr.bonds = substr.bonds!.concat(Array.from(res.bondsToHighLight[j][i]));
        }
      }
      return substr;
    });
  return { rGroups: resCols, highlightCol: substructCol };
}

export async function rGroupsPython(col: DG.Column<string>, core: string, prefix: string, withCore: boolean,
  onlyMatchAtRGroups: boolean):
  Promise<RGroupsRes> {
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
  return { rGroups: resCols };
}
function removeLatestAnalysis(col: DG.Column) {
  if(latestTrellisPlot[col.dataFrame.name] && latestTrellisPlot[col.dataFrame.name]!.dataFrame) 
    latestTrellisPlot[col.dataFrame.name]!.close();
  delete latestTrellisPlot[col.dataFrame.name];
  if(latestAnalysisCols[col.dataFrame.name])
    latestAnalysisCols[col.dataFrame.name]!.forEach((colName: string) => col.dataFrame.columns.remove(colName));
  delete latestAnalysisCols[col.dataFrame.name];
}
