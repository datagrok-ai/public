/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';

import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {addMonomerHoverLink} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';

import {getRules, RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './conversion/pt-rules';
import {doPolyToolConvert} from './conversion/pt-conversion';
import {getOverriddenLibrary} from './conversion/pt-synthetic';
import {defaultErrorHandler} from '../utils/err-info';
import {getLibrariesList} from './utils';
import {
  PT_ERROR_DATAFRAME, PT_UI_ADD_HELM, PT_UI_DIALOG_CONVERSION,
  PT_UI_GET_HELM, PT_UI_LINEARIZE, PT_UI_LINEARIZE_TT,
  PT_UI_HIGHLIGHT_MONOMERS, PT_UI_RULES_USED, PT_UI_HIGHLIGHT_MONOMERS_TT, PT_UI_COLUMN_TT,
} from './const';

import {_package} from '../package';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {MonomerHoverLink} from '@datagrok-libraries/bio/src/monomer-works/utils';
import {ISeqMonomer} from '@datagrok-libraries/bio/src/helm/types';
import {AsmLayoutName, getSeqSrcLayout, MonomerMapProvider} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {addSubstructProvider, getMonomerHover, ISubstruct, setMonomerHover}
  from '@datagrok-libraries/chem-meta/src/types';
import {getMolHighlight} from '@datagrok-libraries/bio/src/monomer-works/seq-to-molfile';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {mergeSubstructs} from '@datagrok-libraries/chem-meta/src/types';
import {dealGroups, helmToMol} from './conversion/pt-atomic';

type PolyToolConvertSerialized = {
  generateHelm: boolean;
  rules: string[];
};

export async function polyToolConvertUI(): Promise<void> {
  await _package.initPromise;
  let dialog: DG.Dialog | null = null;
  try {
    dialog = await getPolyToolConvertDialog();
    dialog?.show();
  } catch (err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.warning('To run PolyTool Conversion, open a dataframe with macromolecules');
    _package.logger.error(errMsg, undefined, errStack);
  }
}

export async function getPolyToolConvertDialog(srcCol?: DG.Column): Promise<DG.Dialog | null> {
  const subs: Unsubscribable[] = [];
  const destroy = () => {
    for (const sub of subs) sub.unsubscribe();
  };
  try {
    let srcColVal: DG.Column<string> | undefined = srcCol;
    const srcColList = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
    const customSrcCols = srcColList.filter((col) => {
      const sh = _package.seqHelper.getSeqHandler(col);
      return sh.notation === NOTATION.CUSTOM;
    });
    if (!srcColVal) {
      if (srcColList.length < 1)
        throw new Error(PT_ERROR_DATAFRAME);

      if (customSrcCols.length < 1) {
        const toAtomicLevelFunc = DG.Func.find({package: 'Bio', name: 'toAtomicLevel'})[0];
        if (toAtomicLevelFunc) {
          toAtomicLevelFunc.prepare().edit();
          return null;
        }
        // eslint-disable-next-line max-len
        grok.shell.warning('Polytool requires a macromolecule column with custom notation. \n\nUse Top menu | Bio | Transform | To Atomic Level.');
        return null;
      }

      srcColVal = srcColList[0];
    }
    const srcColInput = ui.input.column('Column', {
      table: srcColVal.dataFrame, value: srcColVal,
      filter: (col: DG.Column) => {
        return customSrcCols.includes(col);
      }
    });

    const generateHelmInput = ui.input.bool(PT_UI_GET_HELM, {value: true});
    ui.tooltip.bind(generateHelmInput.root, PT_UI_ADD_HELM);

    const linearizeInput = ui.input.bool(PT_UI_LINEARIZE, {value: true});
    ui.tooltip.bind(linearizeInput.root, PT_UI_LINEARIZE_TT);

    const highlightMonomersInput = ui.input.bool(PT_UI_HIGHLIGHT_MONOMERS, {value: false});
    ui.tooltip.bind(highlightMonomersInput.root, PT_UI_HIGHLIGHT_MONOMERS_TT);
    ui.tooltip.bind(srcColInput.root, PT_UI_COLUMN_TT);
    let ruleFileList: string[];
    const ruleInputs = new RuleInputs(RULES_PATH, RULES_STORAGE_NAME, '.json', {
      onValueChanged: (value: string[]) => { ruleFileList = value; }
    });
    const rulesHeader = ui.inlineText([PT_UI_RULES_USED]);
    ui.tooltip.bind(rulesHeader, 'Add or specify rules to use');
    const rulesForm = await ruleInputs.getForm();

    const div = ui.divV([
      srcColInput,
      generateHelmInput,
      linearizeInput,
      highlightMonomersInput,
      rulesHeader,
      rulesForm
    ]);

    const exec = async (): Promise<void> => {
      try {
        const ruleFileList = await ruleInputs.getActive();
        // chirality: no input any more, stereo is always absolute
        await polyToolConvert(
          srcColInput.value!, generateHelmInput.value!, linearizeInput.value!,
          true, highlightMonomersInput.value!, ruleFileList,
        );
      } catch (err: any) {
        defaultErrorHandler(err);
      }
    };

    const dialog = ui.dialog(PT_UI_DIALOG_CONVERSION)
      .add(div)
      .onOK(() => { exec(); });
    subs.push(dialog.onClose.subscribe(() => {
      destroy();
    }));
    dialog.history(
      /* getInput */ (): PolyToolConvertSerialized => {
        return {
          generateHelm: generateHelmInput.value,
          rules: ruleFileList,
        };
      },
      /* applyInput */ (x: PolyToolConvertSerialized): void => {
        generateHelmInput.value = x.generateHelm;
        ruleInputs.setActive(x.rules);
      });
    return dialog;
  } catch (err: any) {
    destroy(); // on failing to build a dialog
    throw err;
  }
}

/** Returns Helm and molfile columns.  */
export async function polyToolConvert(seqCol: DG.Column<string>,
  generateHelm: boolean, linearize: boolean, chiralityEngine: boolean, highlight: boolean, ruleFiles: string[]
): Promise<[DG.Column, DG.Column]> {
  const pi = DG.TaskBarProgressIndicator.create('PolyTool converting...');
  try {
    const getUnusedName = (df: DG.DataFrame | undefined, colName: string): string => {
      if (!df) return colName;
      return df.columns.getUnusedName(colName);
    };
    const helmHelper = await getHelmHelper(); // initializes JSDraw and org

    const table = seqCol.dataFrame;
    const rules = await getRules(ruleFiles);
    const [resList, isLinear, positionMaps] = doPolyToolConvert(seqCol.toList(), rules, helmHelper);

    const resHelmColName = getUnusedName(table, `transformed(${seqCol.name})`);
    const resHelmCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, resHelmColName, resList.length)
      .init((rowIdx: number) => { return resList[rowIdx]; });
    resHelmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    resHelmCol.meta.units = NOTATION.HELM;
    resHelmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');
    if (generateHelm && table) table.columns.add(resHelmCol, true);


    const seqHelper: ISeqHelper = await getSeqHelper();

    const lib = await getOverriddenLibrary(rules);
    const resHelmColTemp = resHelmCol.temp;
    resHelmColTemp[MmcrTemps.overriddenLibrary] = lib;
    resHelmCol.temp = resHelmColTemp;

    const resMolCol = await helmToMol(resHelmCol, resList,
      isLinear, chiralityEngine, highlight, linearize, lib, null, seqHelper);
    resMolCol.name = getUnusedName(table, `molfile(${seqCol.name})`);
    resMolCol.semType = DG.SEMTYPE.MOLECULE;

    if (table) {
      table.columns.add(resMolCol, true);
      await grok.data.detectSemanticTypes(table);
    }

    //buildMonomerHoverLink(resHelmCol, resMolCol, lib, seqHelper, null);
    // Rows redone with the linear layout (linearize) get linear maps, the others non-linear ones
    const rowLayout = (rowIdx: number): AsmLayoutName => linearize && isLinear[rowIdx] ?
      'linear' : getSeqSrcLayout(resMolCol, 'nonlinear');
    buildCyclizedMonomerHoverLink(seqCol, resHelmCol, resMolCol, lib, seqHelper, rowLayout, positionMaps);

    return [resHelmCol, resMolCol];
  } finally {
    pi.close();
  }
}

/** Hover from the cyclized notation column to the molecules made from its HELM: the maps are the assembler's,
 * per HELM row (MonomerMapProvider), and {@link positionMaps} takes a cyclized position to its HELM positions.
 * @param {Function} rowLayout the layout each row's molecule was assembled with */
function buildCyclizedMonomerHoverLink(
  cyclizedCol: DG.Column<string>, seqCol: DG.Column<string>, molCol: DG.Column<string>,
  monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, rowLayout: (rowIdx: number) => AsmLayoutName,
  positionMaps: number[][][]
): MonomerHoverLink {
  const maps = new MonomerMapProvider(seqCol, monomerLib, seqHelper, rowLayout);
  let lastGrid: DG.Grid | null = null;
  let invalidateTimer: number | null = null;
  /** The whole-molecule highlight's maps arrive one by one; one redraw for those arriving together */
  const scheduleInvalidate = (): void => {
    if (invalidateTimer !== null) return;
    invalidateTimer = window.setTimeout(() => {
      invalidateTimer = null;
      lastGrid?.invalidate();
      for (const tv of grok.shell.tableViews)
        if (molCol.dataFrame && tv.dataFrame?.id === molCol.dataFrame.id) tv.grid?.invalidate();
    }, 50);
  };

  const resLink: MonomerHoverLink = {
    targetCol: molCol,
    handler: (seqGridCell: DG.GridCell, cyclizedMonomer: ISeqMonomer | null, targetGridCol: DG.GridColumn): boolean => {
      if (!seqGridCell || !targetGridCol.grid || !seqCol.dataFrame)
        return true;
      const grid = targetGridCol.grid;
      lastGrid = grid;
      const tableRowIdx = seqGridCell.tableRowIndex!;
      const gridRowIdx = seqGridCell.gridRow;
      const targetGridCell = grid.cell(targetGridCol.name, gridRowIdx);
      const positionMap = positionMaps[gridRowIdx];

      const prev = getMonomerHover();
      if (!prev || (prev && (prev.dataFrameId != seqCol.dataFrame?.id || prev.gridRowIdx != gridRowIdx ||
        prev.seqColName != seqCol.name || prev.seqPosition != cyclizedMonomer?.position))
      ) {
        if (prev) {
          setMonomerHover(null);
          //prev.gridCell.grid?.invalidate();
          prev.gridCell.render();
        }
        if (!cyclizedMonomer) {
          setMonomerHover(null);
          return true;
        }

        setMonomerHover({
          gridCell: targetGridCell,
          dataFrameId: seqCol.dataFrame.id,
          gridRowIdx: gridRowIdx,
          seqColName: seqCol.name,
          seqPosition: cyclizedMonomer ? cyclizedMonomer.position : -1,
          getSubstruct: (): ISubstruct | undefined => { // Gets monomer highlight
            if (!cyclizedMonomer || cyclizedMonomer.symbol === '*')
              return undefined;

            const molMonomerMap = maps.get(tableRowIdx); // undefined until the map arrives
            if (!molMonomerMap)
              return undefined;

            const resSubstructList: ISubstruct[] = [];
            const seqMonomerList: number[] = positionMap[cyclizedMonomer.position];
            for (const seqMonomer of seqMonomerList) {
              const monomerMap = molMonomerMap.get(seqMonomer); // single monomer
              if (!monomerMap) return {atoms: [], bonds: [], highlightAtomColors: [], highlightBondColors: []};
              resSubstructList.push(getMolHighlight([monomerMap], monomerLib));
            }
            //TODO: refine merge substract
            const res: ISubstruct = mergeSubstructs(resSubstructList);
            return res;
          }
        });

        // Request the hovered row's map first; redraw the molecule cell when it arrives, if still hovered
        const hoveredPosition = cyclizedMonomer.position;
        maps.get(tableRowIdx, () => {
          const cur = getMonomerHover();
          if (cur && cur.dataFrameId == seqCol.dataFrame?.id && cur.gridRowIdx == gridRowIdx &&
            cur.seqColName == seqCol.name && cur.seqPosition == hoveredPosition)
            targetGridCell.render();
        }, true);

        // TODO: Invalidate targetGridCell
        //grid.invalidate();
        targetGridCell.render();
      }

      return true;
    },
    /* ISubstructProvider.*/getSubstruct: (tableRowIdx: number | null): ISubstruct | undefined =>{
      // Gets whole molecule highlight
      if (molCol.getTag(ChemTags.SEQUENCE_SRC_HL_MONOMERS) != 'true') return undefined;
      if (tableRowIdx == null) return undefined;
      const seq = seqCol.get(tableRowIdx);
      if (!seq) return undefined;

      // Prefetch: the rows the renderer draws get their maps requested; their grids redraw when they arrive
      const molMonomerMap = maps.get(tableRowIdx, scheduleInvalidate);
      if (!molMonomerMap) return undefined;
      const res: ISubstruct = getMolHighlight(molMonomerMap.values(), monomerLib);
      return res;
    }
  };

  addMonomerHoverLink(cyclizedCol.temp, resLink);
  addSubstructProvider(molCol.temp, resLink);

  return resLink;
}
