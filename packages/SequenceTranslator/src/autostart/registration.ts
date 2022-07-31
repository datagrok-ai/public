import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {siRnaAxolabsToGcrs, gcrsToNucleotides, asoGapmersBioSpringToGcrs, gcrsToMermade12,
} from '../structures-works/converters';
import {map, COL_NAMES, MODIFICATIONS} from '../structures-works/map';
import {getFormat} from '../structures-works/sequence-codes-tools';
import {sequenceToMolV3000} from '../structures-works/from-monomers';

import {SALTS_CSV} from '../salts';
import {USERS_CSV} from '../users';
import {ICDS} from '../ICDs';
import {SOURCES} from '../sources';
import {IDPS} from '../IDPs';

const weightsObj: {[code: string]: number} = {};
for (const synthesizer of Object.keys(map)) {
  for (const technology of Object.keys(map[synthesizer])) {
    for (const code of Object.keys(map[synthesizer][technology]))
      weightsObj[code] ?? map[synthesizer][technology][code].weight;
  }
}
for (const [key, value] of Object.entries(MODIFICATIONS))
  weightsObj[key] = value.molecularWeight;


function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}

function stringify(items: string[]): string {
  return '["' + items.join('", "') + '"]';
}

function molecularWeight(sequence: string, weightsObj: {[index: string]: number}): number {
  const codes = sortByStringLengthInDescendingOrder(Object.keys(weightsObj)).concat(Object.keys(MODIFICATIONS));
  let weight = 0;
  let i = 0;
  while (i < sequence.length) {
    const matchedCode = codes.find((s) => s == sequence.slice(i, i + s.length))!;
    weight += weightsObj[sequence.slice(i, i + matchedCode.length)];
    i += matchedCode!.length;
  }
  return weight - 61.97;
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
    const format = getFormat(structureColumn.get(i));
    result += (typeColumn.get(i) == 'SS') ?
      sequenceToMolV3000(structureColumn.get(i), false, true, format!) + '\n' + `>  <Sequence>\nSense Strand\n\n` :
      sequenceToMolV3000(structureColumn.get(i), true, true, format!) + '\n' + `>  <Sequence>\nAnti Sense\n\n`;
    for (const col of table.columns) {
      if (col.name != COL_NAMES.SEQUENCE)
        result += `>  <${col.name}>\n${col.get(i)}\n\n`;
    }
    result += '$$$$\n\n';
  }
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}

//tags: autostart
export function autostartOligoSdFileSubscription() {
  grok.events.onViewAdded.subscribe((v: any) => {
    if (v.type == 'TableView') {
      if (v.dataFrame.columns.contains(COL_NAMES.TYPE))
        oligoSdFile(v.dataFrame);
      grok.events.onContextMenu.subscribe((args) => {
        const seqCol = args.args.context.table.currentCol;
        if (DG.Detector.sampleCategories(seqCol, (s) => /^[fsACGUacgu]{6,}$/.test(s))) {
          args.args.menu.item('Convert Axolabs to GCRS', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaAxolabsToGcrs(seqCol.get(i));
            });
          });
        } else if (DG.Detector.sampleCategories(seqCol, (s) => /^[fmpsACGU]{6,}$/.test(s)) ||
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
          });
        } else if (DG.Detector.sampleCategories(seqCol, (s) => /^[*56789ATGC]{6,}$/.test(s))) {
          args.args.menu.item('Convert Biospring to GCRS', () => {
            const seqCol = args.args.context.table.currentCol;
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return asoGapmersBioSpringToGcrs(seqCol.get(i));
            });
          });
        } else if (DG.Detector.sampleCategories(seqCol, (s) => /^[*1-8]{6,}$/.test(s))) {
          args.args.menu.item('Convert Biospring to GCRS', () => {
            args.args.context.table.columns.addNewString(seqCol.name + ' to GCRS').init((i: number) => {
              return siRnaAxolabsToGcrs(seqCol.get(i));
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

  function addColumns(t: DG.DataFrame, saltsDf: DG.DataFrame) {
    if (t.columns.contains(COL_NAMES.COMPOUND_NAME))
      return grok.shell.error('Columns already exist!');

    const sequence = t.getCol(COL_NAMES.SEQUENCE);
    const salt = t.getCol(COL_NAMES.SALT);
    const equivalents = t.getCol(COL_NAMES.EQUIVALENTS);

    t.columns.addNewString(COL_NAMES.COMPOUND_NAME).init((i: number) => sequence.get(i));
    t.columns.addNewString(COL_NAMES.COMPOUND_COMMENTS).init((i: number) => (i > 0 && i % 2 == 0) ?
      sequence.getString(i) + '; duplex of SS: ' + sequence.getString(i - 2) + ' and AS: ' + sequence.getString(i - 1) :
      sequence.getString(i),
    );
    const molWeightCol = saltsDf.getCol('MOLWEIGHT');
    const saltNamesList = saltsDf.getCol('DISPLAY').toList();
    t.columns.addNewFloat(COL_NAMES.CPD_MW)
      .init((i: number) => molecularWeight(sequence.get(i), weightsObj));
    t.columns.addNewFloat(COL_NAMES.SALT_MASS).init((i: number) => {
      const saltRowIndex = saltNamesList.indexOf(salt.get(i));
      const mw = molWeightCol.get(saltRowIndex);
      return mw * equivalents.get(i);
    });
    t.columns.addNewCalculated(COL_NAMES.BATCH_MW,
      '${' + COL_NAMES.CPD_MW + '} + ${' + COL_NAMES.SALT_MASS + '}', DG.COLUMN_TYPE.FLOAT, false,
    );

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
          addColumns(table, saltsDf);
          grok.shell.tableView(table.name).grid.columns.setOrder(Object.values(COL_NAMES));
        }, 'Add columns: \'' + [COL_NAMES.COMPOUND_NAME, COL_NAMES.COMPOUND_COMMENTS, COL_NAMES.CPD_MW,
          COL_NAMES.SALT_MASS, COL_NAMES.BATCH_MW].join('\', \''), ''),
        ui.button('Save SD file', () => saveTableAsSdFile(addColumnsPressed ? newDf : table)),
      );

      const view = grok.shell.getTableView(table.name);

      view.dataFrame.getCol(COL_NAMES.TYPE).setTag(DG.TAGS.CHOICES, '["AS", "SS", "Duplex"]');
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
          });
        }
      });
    }),
  ]);
  grok.shell.v.setRibbonPanels([[d]]);
}
