import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  siRnaBioSpringToGcrs, siRnaAxolabsToGcrs, gcrsToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12,
  siRnaNucleotidesToGcrs
} from '../hardcode-to-be-eliminated/converters';
import {weightsObj, SYNTHESIZERS} from '../hardcode-to-be-eliminated/map';
import {SEQUENCE_TYPES, COL_NAMES, GENERATED_COL_NAMES} from './constants';
import {saltMass, saltMolWeigth, molecularWeight, batchMolWeight} from './calculations';
import {isValidSequence} from '../sdf-tab/sequence-codes-tools';
import {sequenceToMolV3000} from '../utils/structures-works/from-monomers';
import {linkStrandsV3000} from '../utils/structures-works/mol-transformations';
import {stringify, download, removeEmptyRows, differenceOfTwoArrays} from '../utils/helpers';

import {SALTS_CSV} from '../hardcode-to-be-eliminated/salts';
import {USERS_CSV} from '../hardcode-to-be-eliminated/users';
import {ICDS} from '../hardcode-to-be-eliminated/ICDs';
import {SOURCES} from '../hardcode-to-be-eliminated/sources';
import {IDPS} from '../hardcode-to-be-eliminated/IDPs';

import {sdfAddColumns} from '../utils/sdf-add-columns';
import {sdfSaveTable} from '../utils/sdf-save-table';

const enum PREFIXES {
  AS = 'AS',
  SS = 'SS',
  AS1 = 'AS1',
  AS2 = 'AS2'
}

const enum SEQ_TYPE {
  AS = 'AS',
  SS = 'SS',
  DUPLEX = 'Duplex',
  DIMER = 'Dimer',
}

/** Computable categories of sequence types */
const enum SEQ_TYPE_CATEGORY {
  AS_OR_SS,
  DUPLEX,
  DIMER,
}

/** Map between types and their categories inferrable from 'Sequence' column */
const typeCategoryMap = {
  [SEQ_TYPE.AS]: SEQ_TYPE_CATEGORY.AS_OR_SS,
  [SEQ_TYPE.SS]: SEQ_TYPE_CATEGORY.AS_OR_SS,
  [SEQ_TYPE.DIMER]: SEQ_TYPE_CATEGORY.DIMER,
  [SEQ_TYPE.DUPLEX]: SEQ_TYPE_CATEGORY.DUPLEX,
};

/** Style used for cells in 'Type' column  */
const typeColCellStyle = {
  'display': 'flex',
  'justify-content': 'center',
  'align-items': 'center',
  'text-color': 'var(--grey-5)', // --grey-6 does not match other cells
  'width': '100%',
  'height': '100%',
};

const pinkBackground = {
  'background-color': '#ff8080',
};

/** Style used for a cell with invalid value  */
const typeColErrorStyle = Object.assign({}, pinkBackground, typeColCellStyle);

export function sdfHandleErrorUI(msgPrefix: string, df: DG.DataFrame, rowI: number, err: any) {
  const errStr: string = err.toString();
  const errMsg: string = msgPrefix + `row #${rowI + 1}, name: '${df.get('Chemistry Name', rowI)}', ` +
    `type: ${df.get('Type', rowI)} error: ${errStr}.`;
  grok.shell.warning(errMsg);
}

/** Determine the category of the value specified in 'Types' column  */
function getActualTypeClass(actualType: string): SEQ_TYPE_CATEGORY {
  if (Object.keys(typeCategoryMap).includes(actualType))
    return typeCategoryMap[actualType as SEQ_TYPE];
  else
    throw new Error('Some types in \'Types\' column are invalid ');
}

function isASorSS(splittedLines: string[][]): boolean {
  return splittedLines.length === 1 && splittedLines[0].length === 1;
}

/** Check whether the number of lines and prefixes in the 'Sequence' string
 * are valid  */
function verifyPrefixes(splittedLines: string[][], allowedPrefixes: Set<PREFIXES>, allowedLength: number): boolean {
  const lengthCriterion = splittedLines.length === allowedLength;
  let prefixCriterion = true;
  for (const line of splittedLines) {
    const prefix = line[0];
    prefixCriterion &&= (allowedPrefixes.has(prefix as PREFIXES));
  }
  return lengthCriterion && prefixCriterion;
}

function isDuplex(splittedLines: string[][]): boolean {
  const allowedPrefixes = new Set([PREFIXES.SS, PREFIXES.AS]);
  return verifyPrefixes(splittedLines, allowedPrefixes, 2);
}

function isDimer(splittedLines: string[][]): boolean {
  const allowedPrefixes = new Set([PREFIXES.SS, PREFIXES.AS1, PREFIXES.AS2]);
  return verifyPrefixes(splittedLines, allowedPrefixes, 3);
}

function inferTypeClassFromSequence(seq: string): SEQ_TYPE_CATEGORY {
  const lines = seq.split('\n');
  const splittedLines = [];
  for (const line of lines)
    splittedLines.push(line.split(' '));
  if (isASorSS(splittedLines))
    return SEQ_TYPE_CATEGORY.AS_OR_SS;
  else if (isDuplex(splittedLines))
    return SEQ_TYPE_CATEGORY.DUPLEX;
  else if (isDimer(splittedLines))
    return SEQ_TYPE_CATEGORY.DIMER;
  else
    throw new Error('Some cells in \'Sequence\' column have wrong formatting');
}

/** Compare type specified in 'Type' column to that computed from 'Sequence' column  */
function validateType(actualType: string, seq: string): boolean {
  if (actualType === '' && seq === '')
    return true;
  else
    return getActualTypeClass(actualType) === inferTypeClassFromSequence(seq);
}

function oligoSdFileGrid(view: DG.TableView): void {
  const typeColName = 'Type';
  const seqColName = 'Sequence';
  const grid = view.grid;
  const df = view.dataFrame;
  const typeCol = df.getCol(typeColName);
  grid.columns.byName(typeColName)!.cellType = 'html';
  const seqCol = df.getCol(seqColName);
  grid.onCellPrepare((gridCell: DG.GridCell) => {
    if (gridCell.isTableCell && gridCell.gridColumn.column!.name === typeColName) {
      let isValidType = false;
      let formattingError = false;
      try {
        isValidType = validateType(gridCell.cell.value, seqCol.get(gridCell.tableRow!.idx));
      } catch {
        formattingError = true;
      }
      const el = ui.div(
        gridCell.cell.value, isValidType ? {style: typeColCellStyle} : {style: typeColErrorStyle}
      );
      gridCell.style.element = el;
      const msg = formattingError ? 'Sequence pattern or Type value has wrong formatting' :
        'Input in Type column doesn\'t match the Sequence pattern';
      if (!isValidType)
        ui.tooltip.bind(el, msg);
    }
  });
}

export function autostartOligoSdFileSubscription() {
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type === DG.VIEW_TYPE.TABLE_VIEW) {
      if (v.dataFrame.columns.contains(COL_NAMES.TYPE)) {
        oligoSdFileGrid(v);
        oligoSdFile(v.dataFrame);
      }

      // Should be removed after fixing bug https://github.com/datagrok-ai/public/issues/808
      grok.events.onContextMenu.subscribe((args) => {
        if (!(args.args.context instanceof DG.Grid)) return;
        const grid: DG.Grid = args.args.context as DG.Grid;
        const menu: DG.Menu = args.args.menu;

        const seqCol = grid.table.currentCol; // /^[fsACGUacgu]{6,}$/
        if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C){6,}$/.test(s))) {
          menu.item('Convert raw nucleotides to GCRS', () => {
            grid.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaNucleotidesToGcrs(seqCol.get(i));
            });
          });
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|s|A|C|G|U|a|c|g|u){6,}$/.test(s))) {
          menu.item('Convert Axolabs to GCRS', () => {
            grid.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaAxolabsToGcrs(seqCol.get(i));
            });
          }); // /^[fmpsACGU]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|m|ps|A|C|G|U){6,}$/.test(s)) ||
          DG.Detector.sampleCategories(seqCol, (s) => /^(?=.*moe)(?=.*5mC)(?=.*ps){6,}/.test(s))) {
          menu.item('Convert GCRS to raw', () => {
            grid.table.columns.addNewString(seqCol.name + ' to raw').init((i: number) => {
              return gcrsToNucleotides(seqCol.get(i));
            });
          });
          menu.item('Convert GCRS to MM12', () => {
            grid.table.columns.addNewString(seqCol.name + ' to MM12').init((i: number) => {
              return gcrsToMermade12(seqCol.get(i));
            });
          }); // /^[*56789ATGC]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|5|6|7|8|9|A|T|G|C){6,}$/.test(s))) {
          menu.item('Convert Biospring to GCRS', () => {
            const seqCol = grid.table.currentCol;
            grid.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return asoGapmersBioSpringToGcrs(seqCol.get(i));
            });
          }); // /^[*1-8]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|1|2|3|4|5|6|7|8){6,}$/.test(s))) {
          menu.item('Convert Biospring to GCRS', () => {
            grid.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaBioSpringToGcrs(seqCol.get(i));
            });
          });
        }
      });
    }
  });
}

export function oligoSdFile(table: DG.DataFrame) {
  const saltsDf = DG.DataFrame.fromCsv(SALTS_CSV);
  const usersDf = DG.DataFrame.fromCsv(USERS_CSV);
  const sourcesDf = DG.DataFrame.fromCsv(SOURCES);
  const icdsDf = DG.DataFrame.fromCsv(ICDS);
  const idpsDf = DG.DataFrame.fromCsv(IDPS);

  const saltCol = table.getCol(COL_NAMES.SALT);
  const equivalentsCol = table.getCol(COL_NAMES.EQUIVALENTS);

  const saltsMolWeightList: number[] = saltsDf.getCol('MOLWEIGHT').toList();
  const saltNamesList: string[] = saltsDf.getCol('DISPLAY').toList();

  let newDf: DG.DataFrame | undefined = undefined;

  const d = ui.div([
    ui.icons.edit(() => {
      d.innerHTML = '';
      if (table.getCol(COL_NAMES.IDP).type != DG.COLUMN_TYPE.STRING)
        table.changeColumnType(COL_NAMES.IDP, DG.COLUMN_TYPE.STRING);
      d.append(
        ui.divH([
          ui.button('Add columns',
            () => {
              newDf = sdfAddColumns(table, saltNamesList, saltsMolWeightList,
                (rowI, err) => { sdfHandleErrorUI('Error on ', table, rowI, err); });
              grok.shell.getTableView(newDf.name).grid.columns.setOrder(Object.values(COL_NAMES));
            },
            `Add columns: '${GENERATED_COL_NAMES.join(`', '`)}'`),
          ui.bigButton('Save SDF', () => {
            const df: DG.DataFrame = newDf ?? table;
            sdfSaveTable(df,
              (rowI, err) => { sdfHandleErrorUI('Skip ', df, rowI, err); });
          }, 'Save SD file'),
        ])
      );

      const view = grok.shell.getTableView(table.name);
      view.grid.setOptions({rowHeight: 45});
      view.dataFrame.getCol(COL_NAMES.TYPE).setTag(DG.TAGS.CHOICES, stringify(Object.values(SEQUENCE_TYPES)));
      view.dataFrame.getCol(COL_NAMES.OWNER).setTag(DG.TAGS.CHOICES, stringify(usersDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.SALT).setTag(DG.TAGS.CHOICES, stringify(saltsDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.SOURCE).setTag(DG.TAGS.CHOICES, stringify(sourcesDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.ICD).setTag(DG.TAGS.CHOICES, stringify(icdsDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.IDP).setTag(DG.TAGS.CHOICES, stringify(idpsDf.columns.byIndex(0).toList()));

      grok.events.onContextMenu.subscribe((args) => {
        if (!(args.args.context instanceof DG.Grid)) return;
        const grid: DG.Grid = args.args.context as DG.Grid;
        const menu: DG.Menu = args.args.menu;

        if ([COL_NAMES.TYPE, COL_NAMES.OWNER, COL_NAMES.SALT, COL_NAMES.SOURCE, COL_NAMES.ICD, COL_NAMES.IDP]
          .includes(grid.table.currentCol.name)) {
          menu.item('Fill Column With Value', () => {
            const v = grid.table.currentCell.value;
            grid.table.currentCell.column.init(v);
            for (let i = 0; i < view.dataFrame.rowCount; i++)
              updateCalculatedColumns(view.dataFrame, i);
          });
        }
      });

      view.dataFrame.onDataChanged.subscribe(() => {
        const colName = view.dataFrame.currentCol.name;
        if ([COL_NAMES.SALT, COL_NAMES.EQUIVALENTS, COL_NAMES.SALT_MOL_WEIGHT].includes(colName))
          updateCalculatedColumns(view.dataFrame, view.dataFrame.currentRowIdx);
      });

      function updateCalculatedColumns(t: DG.DataFrame, i: number): void {
        const smValue = saltMass(saltNamesList, saltsMolWeightList, equivalentsCol, i, saltCol);
        t.getCol(COL_NAMES.SALT_MASS).set(i, smValue, false);
        const smwValue = saltMolWeigth(saltNamesList, saltCol, saltsMolWeightList, i);
        t.getCol(COL_NAMES.SALT_MOL_WEIGHT).set(i, smwValue, false);
        const bmw = batchMolWeight(t.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT), t.getCol(COL_NAMES.SALT_MASS), i);
        t.getCol(COL_NAMES.BATCH_MOL_WEIGHT).set(i, bmw, false);
      }
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
