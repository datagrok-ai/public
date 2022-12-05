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

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {parseStrandsFromDuplexCell, parseStrandsFromTriplexOrDimerCell} from '../utils/parse';

async function saveTableAsSdFile(table: DG.DataFrame) {
  if (GENERATED_COL_NAMES.some((colName) => !table.columns.contains(colName))) {
    const absentColNames = differenceOfTwoArrays(GENERATED_COL_NAMES, table.columns.names()).join(`', '`);
    grok.shell.warning(`File saved without columns '${absentColNames}'`);
  }

  const sequenceCol = table.getCol(COL_NAMES.SEQUENCE);
  const typeCol = table.getCol(COL_NAMES.TYPE);

  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    const format = SYNTHESIZERS.GCRS; //getFormat(sequenceCol.get(i))!;
    if (typeCol.get(i) == SEQUENCE_TYPES.SENSE_STRAND)
      result += `${sequenceToMolV3000(sequenceCol.get(i), false, true, format)}\n> <Sequence>\nSense Strand\n\n`;
    else if (typeCol.get(i) == SEQUENCE_TYPES.ANTISENSE_STRAND)
      result += `${sequenceToMolV3000(sequenceCol.get(i), true, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
    else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
      const obj = parseStrandsFromDuplexCell(sequenceCol.get(i));
      const as = `${sequenceToMolV3000(obj.AS, true, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
      const ss = `${sequenceToMolV3000(obj.SS, false, true, format)}\n> <Sequence>\nSense Strand\n\n`;
      result += `${linkStrandsV3000({senseStrands: [ss], antiStrands: [as]}, true)}\n\n`;
    } else if ([SEQUENCE_TYPES.TRIPLEX, SEQUENCE_TYPES.DIMER].includes(typeCol.get(i))) {
      const obj = parseStrandsFromTriplexOrDimerCell(sequenceCol.get(i));
      const as1 = `${sequenceToMolV3000(obj.AS1, true, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
      const as2 = `${sequenceToMolV3000(obj.AS2, true, true, format)}\n> <Sequence>\nAnti Sense\n\n`;
      const ss = `${sequenceToMolV3000(obj.SS, false, true, format)}\n> <Sequence>\nSense Strand\n\n`;
      result += `${linkStrandsV3000({senseStrands: [ss], antiStrands: [as1, as2]}, true)}\n\n`;
    }

    for (const col of table.columns) {
      if (col.name != COL_NAMES.SEQUENCE)
        result += `> <${col.name}>\n${col.get(i)}\n\n`;
    }
    result += '$$$$\n';
  }
  download(`${table.name}.sdf`, encodeURIComponent(result));
}

export function autostartOligoSdFileSubscription() {
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type == DG.VIEW_TYPE.TABLE_VIEW) {
      if (v.dataFrame.columns.contains(COL_NAMES.TYPE))
        oligoSdFile(v.dataFrame);

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

  const sequenceCol = table.getCol(COL_NAMES.SEQUENCE);
  const saltCol = table.getCol(COL_NAMES.SALT);
  const equivalentsCol = table.getCol(COL_NAMES.EQUIVALENTS);
  const typeCol = table.getCol(COL_NAMES.TYPE);
  const chemistryNameCol = table.getCol(COL_NAMES.CHEMISTRY_NAME);

  const molWeightCol = saltsDf.getCol('MOLWEIGHT');
  const saltNamesList = saltsDf.getCol('DISPLAY').toList();

  let newDf: DG.DataFrame;
  let addColumnsPressed = false;

  function addColumns(t: DG.DataFrame) {
    if (GENERATED_COL_NAMES.some((colName) => t.columns.contains(colName)))
      return grok.shell.error('Columns already exist');

    t = removeEmptyRows(t, sequenceCol);

    function handleError(rowI: number, err: any) {
      const errStr: string = err.toString();
      const errMsg: string = `Table row #${rowI + 1}, name: '${t.get('Chemistry Name', rowI)}', ` +
        `type: ${t.get('Type', rowI)} error: ${errStr}.`;
      grok.shell.warning(errMsg);
    }

    t.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => {
      let res: string = '';
      try {
        res = ([SEQUENCE_TYPES.DUPLEX, SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) ?
          chemistryNameCol.get(i) :
          sequenceCol.get(i);
      } catch (err) {
        handleError(i, err);
      }
      return res;
    });

    t.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => {
      let res: string = '';
      try {
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i)))
          res = sequenceCol.get(i);
        else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
          const obj = parseStrandsFromDuplexCell(sequenceCol.get(i));
          res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.SS} and AS: ${obj.AS}`;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
          const obj = parseStrandsFromTriplexOrDimerCell(sequenceCol.get(i));
          res = `${chemistryNameCol.get(i)}; duplex of SS: ${obj.SS} and AS1: ${obj.AS1} and AS2: ${obj.AS2}`;
        }
      } catch (err) {
        handleError(i, err);
      }
      return res;
    });

    t.columns.addNewFloat(COL_NAMES.COMPOUND_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        if ([SEQUENCE_TYPES.SENSE_STRAND, SEQUENCE_TYPES.ANTISENSE_STRAND].includes(typeCol.get(i))) {
          res = (isValidSequence(sequenceCol.get(i), null).indexOfFirstNotValidChar == -1) ?
            molecularWeight(sequenceCol.get(i), weightsObj) :
            DG.FLOAT_NULL;
        } else if (typeCol.get(i) == SEQUENCE_TYPES.DUPLEX) {
          const obj = parseStrandsFromDuplexCell(sequenceCol.get(i));
          res = (Object.values(obj).every((seq) => isValidSequence(seq, null).indexOfFirstNotValidChar == -1)) ?
            molecularWeight(obj.SS, weightsObj) + molecularWeight(obj.AS, weightsObj) :
            DG.FLOAT_NULL;
        } else if ([SEQUENCE_TYPES.DIMER, SEQUENCE_TYPES.TRIPLEX].includes(typeCol.get(i))) {
          const obj = parseStrandsFromTriplexOrDimerCell(sequenceCol.get(i));
          res = (Object.values(obj).every((seq) => isValidSequence(seq, null).indexOfFirstNotValidChar == -1)) ?
            molecularWeight(obj.SS, weightsObj) + molecularWeight(obj.AS1, weightsObj) +
            molecularWeight(obj.AS2, weightsObj) :
            DG.FLOAT_NULL;
        }
      } catch (err) {
        handleError(i, err);
      }
      return res;
    });

    t.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = saltMass(saltNamesList, molWeightCol, equivalentsCol, i, saltCol);
      } catch (err) {
        handleError(i, err);
      }
      return res;
    });

    t.columns.addNewFloat(COL_NAMES.SALT_MOL_WEIGHT).init((i: number) => {
      let res: number = Number.NaN;
      try {
        res = saltMolWeigth(saltNamesList, saltCol, molWeightCol, i);
      } catch (err) {
        handleError(i, err);
      }
      return res;
    });

    const compoundMolWeightCol = t.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT);
    const saltMassCol = t.getCol(COL_NAMES.SALT_MASS);
    t.columns.addNewFloat(COL_NAMES.BATCH_MOL_WEIGHT).init((i: number) =>
      batchMolWeight(compoundMolWeightCol, saltMassCol, i));

    grok.shell.getTableView(table.name).grid.columns.setOrder(Object.values(COL_NAMES));
    addColumnsPressed = true;
    return newDf = t;
  }

  const d = ui.div([
    ui.icons.edit(() => {
      d.innerHTML = '';
      if (table.getCol(COL_NAMES.IDP).type != DG.COLUMN_TYPE.STRING)
        table.changeColumnType(COL_NAMES.IDP, DG.COLUMN_TYPE.STRING);
      d.append(
        ui.divH([
          ui.button('Add columns', () => addColumns(table), `Add columns: '${GENERATED_COL_NAMES.join(`', '`)}'`),
          ui.bigButton('Save SDF', () => saveTableAsSdFile(addColumnsPressed ? newDf : table), 'Save SD file'),
        ]),
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
        const smValue = saltMass(saltNamesList, molWeightCol, equivalentsCol, i, saltCol);
        t.getCol(COL_NAMES.SALT_MASS).set(i, smValue, false);
        const smwValue = saltMolWeigth(saltNamesList, saltCol, molWeightCol, i);
        t.getCol(COL_NAMES.SALT_MOL_WEIGHT).set(i, smwValue, false);
        const bmw = batchMolWeight(t.getCol(COL_NAMES.COMPOUND_MOL_WEIGHT), t.getCol(COL_NAMES.SALT_MASS), i);
        t.getCol(COL_NAMES.BATCH_MOL_WEIGHT).set(i, bmw, false);
      }
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
