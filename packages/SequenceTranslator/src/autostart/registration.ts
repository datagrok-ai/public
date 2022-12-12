import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  siRnaBioSpringToGcrs, siRnaAxolabsToGcrs, gcrsToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12,
  siRnaNucleotidesToGcrs
} from '../structures-works/converters';
import {weightsObj, SYNTHESIZERS} from '../structures-works/map';
import {SEQUENCE_TYPES, COL_NAMES, GENERATED_COL_NAMES} from './constants';
import {saltMass, saltMolWeigth, molecularWeight, batchMolWeight} from './calculations';
import {isValidSequence} from '../structures-works/sequence-codes-tools';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkStrandsV3000} from '../structures-works/mol-transformations';
import {stringify, download, removeEmptyRows, differenceOfTwoArrays} from '../helpers';

import {SALTS_CSV} from './salts';
import {USERS_CSV} from './users';
import {ICDS} from './ICDs';
import {SOURCES} from './sources';
import {IDPS} from './IDPs';

import {sdfAddColumns} from '../utils/sdf-add-columns';
import {sdfSaveTable} from '../utils/sdf-save-table';

const enum SEQ_TYPE {
  AS = 'AS',
  SS = 'SS',
  DUPLEX = 'Duplex',
  DIMER = 'Dimer',
}

/** Computable classes of sequence types */
const enum SEQ_TYPE_CLASS {
  AS_OR_SS,
  DUPLEX,
  DIMER,
}

/** Style used for a cell with invalid value  */
const errorStyle = {
  'background-color': '#ff8080',
  'width': '100%',
  'height': '100%',
};

export function sdfHandleErrorUI(msgPrefix: string, df: DG.DataFrame, rowI: number, err: any) {
  const errStr: string = err.toString();
  const errMsg: string = msgPrefix + `row #${rowI + 1}, name: '${df.get('Chemistry Name', rowI)}', ` +
    `type: ${df.get('Type', rowI)} error: ${errStr}.`;
  grok.shell.warning(errMsg);
}

// todo: use a dictionary instead?
function getActualTypeClass(actualType: string): SEQ_TYPE_CLASS {
  if (actualType === SEQ_TYPE.AS || actualType === SEQ_TYPE.SS)
    return SEQ_TYPE_CLASS.AS_OR_SS;
  else if (actualType === SEQ_TYPE.DIMER)
    return SEQ_TYPE_CLASS.DIMER;
  else if (actualType === SEQ_TYPE.DUPLEX)
    return SEQ_TYPE_CLASS.DUPLEX;
  else
    throw new Error('Some types in \'Types\' column are invalid ');
}

function inferTypeClassFromSequence(seq: string): SEQ_TYPE_CLASS {
  const lines = seq.split('\n');
  if (lines.length === 1)
    return SEQ_TYPE_CLASS.AS_OR_SS;
  else if (lines.length === 2)
    return SEQ_TYPE_CLASS.DUPLEX;
  else if (lines.length === 3)
    return SEQ_TYPE_CLASS.DIMER;
  else
    throw new Error('Wrong formatting of sequences in \'Sequence\' column');
  //todo: throw in the case of wrong formatting
}

/** Compare type specified in 'Type' column to that computed from 'Sequence' column  */
function validateType(actualType: string, seq: string): boolean {
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
      const isValidType = validateType(gridCell.cell.value, seqCol.get(gridCell.tableRow!.idx));
      gridCell.style.element = ui.div(gridCell.cell.value, isValidType ? {} : {style: errorStyle});
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
