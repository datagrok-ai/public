import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {siRnaBioSpringToGcrs, siRnaAxolabsToGcrs, gcrsToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12,
  siRnaNucleotidesToGcrs} from '../structures-works/converters';
import {map, COL_NAMES, MODIFICATIONS} from '../structures-works/map';
import {isValidSequence} from '../structures-works/sequence-codes-tools';
import {sequenceToMolV3000} from '../structures-works/from-monomers';
import {linkV3000} from '../structures-works/mol-transformations';
import {stringify, sortByStringLengthInDescendingOrder} from '../helpers';

import {SALTS_CSV} from '../salts';
import {USERS_CSV} from '../users';
import {ICDS} from '../ICDs';
import {SOURCES} from '../sources';
import {IDPS} from '../IDPs';

const SEQUENCE_TYPES = {
  SENSE_STRAND: 'SS',
  ANTISENSE_STRAND: 'AS',
  DUPLEX: 'Duplex',
};

function saltMass(saltNames: string[], molWeightCol: DG.Column, equivalentsCol: DG.Column, i: number,
  saltCol: DG.Column) {
  const saltRowIndex = saltNames.indexOf(saltCol.get(i));
  return (
    saltRowIndex == -1 || molWeightCol.get(saltRowIndex) == DG.FLOAT_NULL || equivalentsCol.get(i) == DG.INT_NULL) ?
    DG.FLOAT_NULL :
    molWeightCol.get(saltRowIndex) * equivalentsCol.get(i);
}

function saltMolWeigth(saltNamesList: string[], saltCol: DG.Column, molWeightCol: DG.Column, i: number) {
  const saltRowIndex = saltNamesList.indexOf(saltCol.get(i));
  return (saltRowIndex == -1) ? DG.FLOAT_NULL : molWeightCol.get(saltRowIndex);
}

function batchMolWeight(compoundMolWeightCol: DG.Column, saltMassCol: DG.Column, i: number) {
  return (compoundMolWeightCol.getString(i) == '' || saltMassCol.getString(i) == '') ?
    DG.FLOAT_NULL :
    compoundMolWeightCol.get(i) + saltMassCol.get(i);
}

function molecularWeight(sequence: string, weightsObj: {[index: string]: number}): number {
  const codes = sortByStringLengthInDescendingOrder(Object.keys(weightsObj)).concat(Object.keys(MODIFICATIONS));
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += weightsObj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode.length;
  }
  return weight - 61.97;
}

function parseStrandsFromDuplexCell(s: string): string[] {
  return s.slice(3).split('\r\nAS ');
}

async function saveTableAsSdFile(table: DG.DataFrame) {
  if (!table.columns.contains('Compound Name')) {
    grok.shell.warning(
      'File saved without columns \'' +
      [COL_NAMES.COMPOUND_NAME, COL_NAMES.COMPOUND_COMMENTS, COL_NAMES.CPD_MW,
        COL_NAMES.SALT_MASS, COL_NAMES.BATCH_MW].join('\', \''),
    );
  }
  const structureColumn = table.getCol(COL_NAMES.SEQUENCE);
  const typeColumn = table.getCol(COL_NAMES.TYPE);
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    const format = 'Janssen GCRS Codes'; //getFormat(structureColumn.get(i))!;
    if (typeColumn.get(i) == SEQUENCE_TYPES.DUPLEX) {
      const array = parseStrandsFromDuplexCell(structureColumn.get(i));
      const as = sequenceToMolV3000(array[1], true, true, format) +
      '\n' + `> <Sequence>\nAnti Sense\n\n`;
      const ss = sequenceToMolV3000(array[0], false, true, format) +
      '\n' + `> <Sequence>\nSense Strand\n\n`;
      result += linkV3000([ss, as], true, true) + '\n\n';
    } else if (typeColumn.get(i) == SEQUENCE_TYPES.SENSE_STRAND) {
      const molSS = sequenceToMolV3000(structureColumn.get(i), false, true, format) +
      '\n' + `> <Sequence>\nSense Strand\n\n`;
      result += molSS;
    } else if (typeColumn.get(i) == SEQUENCE_TYPES.ANTISENSE_STRAND) {
      const molAS = sequenceToMolV3000(structureColumn.get(i), true, true, format) +
        '\n' + `> <Sequence>\nAnti Sense\n\n`;
      result += molAS;
    }

    for (const col of table.columns) {
      if (col.name != COL_NAMES.SEQUENCE)
        result += `> <${col.name}>\n${col.get(i)}\n\n`;
    }
    result += '$$$$\n';
  }
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

export function autostartOligoSdFileSubscription() {
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type == 'TableView') {
      if (v.dataFrame.columns.contains(COL_NAMES.TYPE))
        oligoSdFile(v.dataFrame);
      grok.events.onContextMenu.subscribe((args) => {
        const seqCol = args.args.context.table.currentCol; // /^[fsACGUacgu]{6,}$/
        if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|A|U|G|C){6,}$/.test(s))) {
          args.args.menu.item('Convert raw nucleotides to GCRS', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaNucleotidesToGcrs(seqCol.get(i));
            });
          });
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|s|A|C|G|U|a|c|g|u){6,}$/.test(s))) {
          args.args.menu.item('Convert Axolabs to GCRS', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaAxolabsToGcrs(seqCol.get(i));
            });
          }); // /^[fmpsACGU]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|f|m|ps|A|C|G|U){6,}$/.test(s)) ||
        DG.Detector.sampleCategories(seqCol, (s) => /^(?=.*moe)(?=.*5mC)(?=.*ps){6,}/.test(s))) {
          args.args.menu.item('Convert GCRS to raw', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to raw').init((i: number) => {
              return gcrsToNucleotides(seqCol.get(i));
            });
          });
          args.args.menu.item('Convert GCRS to MM12', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to MM12').init((i: number) => {
              return gcrsToMermade12(seqCol.get(i));
            });
          }); // /^[*56789ATGC]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|5|6|7|8|9|A|T|G|C){6,}$/.test(s))) {
          args.args.menu.item('Convert Biospring to GCRS', () => {
            const seqCol = args.args.context.table.currentCol;
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return asoGapmersBioSpringToGcrs(seqCol.get(i));
            });
          }); // /^[*1-8]{6,}$/
        } else if (DG.Detector.sampleCategories(seqCol,
          (s) => /(\(invabasic\)|\(GalNAc-2-JNJ\)|\*|1|2|3|4|5|6|7|8){6,}$/.test(s))) {
          args.args.menu.item('Convert Biospring to GCRS', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
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
  const typeColumn = table.getCol(COL_NAMES.TYPE);
  const chemistryNameCol = table.getCol(COL_NAMES.CHEMISTRY_NAME);

  const molWeightCol = saltsDf.getCol('MOLWEIGHT');
  const saltNamesList = saltsDf.getCol('DISPLAY').toList();

  function addColumns(t: DG.DataFrame) {
    if (t.columns.contains(COL_NAMES.COMPOUND_NAME))
      return grok.shell.error('Columns already exist');

    for (let i = t.rowCount - 1; i > -1; i--) {
      if (sequenceCol.get(i) == '')
        t.rows.removeAt(i, 1, false);
    }

    t.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => {
      return (typeColumn.get(i) == 'Duplex') ? chemistryNameCol.get(i) : sequenceCol.get(i);
    });

    t.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => {
      if (typeColumn.get(i) == 'Duplex') {
        const arr = parseStrandsFromDuplexCell(sequenceCol.get(i));
        return chemistryNameCol.get(i) + '; duplex of SS: ' + arr[0] + ' and AS: ' + arr[1];
      }
      return sequenceCol.get(i);
    });

    const weightsObj: {[code: string]: number} = {};
    for (const synthesizer of Object.keys(map)) {
      for (const technology of Object.keys(map[synthesizer])) {
        for (const code of Object.keys(map[synthesizer][technology]))
          weightsObj[code] = map[synthesizer][technology][code].weight!;
      }
    }
    for (const [key, value] of Object.entries(MODIFICATIONS))
      weightsObj[key] = value.molecularWeight;

    t.columns.addNewFloat(COL_NAMES.CPD_MW).init((i: number) => {
      if (typeColumn.get(i) == 'Duplex') {
        const arr = parseStrandsFromDuplexCell(sequenceCol.get(i));
        return (
          isValidSequence(arr[0], null).indexOfFirstNotValidChar == -1 &&
          isValidSequence(arr[1], null).indexOfFirstNotValidChar == -1
        ) ?
          molecularWeight(arr[0], weightsObj) + molecularWeight(arr[1], weightsObj) :
          DG.FLOAT_NULL;
      }
      return (isValidSequence(sequenceCol.get(i), null).indexOfFirstNotValidChar == -1) ?
        molecularWeight(sequenceCol.get(i), weightsObj) :
        DG.FLOAT_NULL;
    });

    t.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) =>
      saltMass(saltNamesList, molWeightCol, equivalentsCol, i, saltCol));

    t.columns.addNewFloat(COL_NAMES.SALT_MOL_WEIGHT).init((i: number) =>
      saltMolWeigth(saltNamesList, saltCol, molWeightCol, i));

    t.columns.addNewFloat(COL_NAMES.BATCH_MW).init((i: number) =>
      batchMolWeight(t.getCol(COL_NAMES.CPD_MW), t.getCol(COL_NAMES.SALT_MASS), i));

    addColumnsPressed = true;
    return newDf = t;
  }

  let newDf: DG.DataFrame;
  let addColumnsPressed = false;

  const d = ui.div([
    ui.icons.edit(() => {
      d.innerHTML = '';
      if (table.getCol(COL_NAMES.IDP).type != DG.COLUMN_TYPE.STRING)
        table.changeColumnType(COL_NAMES.IDP, DG.COLUMN_TYPE.STRING);
      d.append(
        ui.link('Add Columns', () => {
          addColumns(table);
          view.grid.columns.setOrder(Object.values(COL_NAMES));
        }, 'Add columns: \'' + [COL_NAMES.COMPOUND_NAME, COL_NAMES.COMPOUND_COMMENTS, COL_NAMES.CPD_MW,
          COL_NAMES.SALT_MASS, COL_NAMES.BATCH_MW].join('\', \''), '',
        ),
        ui.button('Save SD file', () => saveTableAsSdFile(addColumnsPressed ? newDf : table)),
      );

      const view = grok.shell.getTableView(table.name);

      view.dataFrame.getCol(COL_NAMES.TYPE).setTag(DG.TAGS.CHOICES, stringify(Object.values(SEQUENCE_TYPES)));
      view.dataFrame.getCol(COL_NAMES.OWNER).setTag(DG.TAGS.CHOICES, stringify(usersDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.SALT).setTag(DG.TAGS.CHOICES, stringify(saltsDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.SOURCE).setTag(DG.TAGS.CHOICES, stringify(sourcesDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.ICD).setTag(DG.TAGS.CHOICES, stringify(icdsDf.columns.byIndex(0).toList()));
      view.dataFrame.getCol(COL_NAMES.IDP).setTag(DG.TAGS.CHOICES, stringify(idpsDf.columns.byIndex(0).toList()));

      grok.events.onContextMenu.subscribe((args) => {
        if ([COL_NAMES.TYPE, COL_NAMES.OWNER, COL_NAMES.SALT, COL_NAMES.SOURCE, COL_NAMES.ICD, COL_NAMES.IDP]
          .includes(args.args.context.table.currentCol.name)) {
          args.args.menu.item('Fill Column With Value', () => {
            const v = args.args.context.table.currentCell.value;
            args.args.context.table.currentCell.column.init(v);
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
        const bmw = batchMolWeight(t.getCol(COL_NAMES.CPD_MW), t.getCol(COL_NAMES.SALT_MASS), i);
        t.getCol(COL_NAMES.BATCH_MW).set(i, bmw, false);
      }
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
