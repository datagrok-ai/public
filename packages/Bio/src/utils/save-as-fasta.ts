import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import wu from 'wu';

import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

import {_package} from '../package';

const FASTA_LINE_WIDTH = 60;

/** Shows dialog to select id columns list and seq column, builds and downloads FASTA content */
export function saveAsFastaUI(): void {
  // Use grid for column order adjusted by user
  const grid: DG.Grid = grok.shell.tv.grid;

  const idGColList: DG.GridColumn[] = wu.count(0).take(grid.columns.length)
    .map((colI: number) => grid.columns.byIndex(colI)!)
    .filter((gcol: DG.GridColumn) => gcol.column ? gcol.column.semType !== DG.SEMTYPE.MACROMOLECULE : false).toArray();
  const defaultIdGCol: DG.GridColumn | undefined = idGColList
    .find((gcol: DG.GridColumn) => gcol.name.toLowerCase().indexOf('id') !== -1);
  const idDefaultValue = defaultIdGCol ? [defaultIdGCol.name] : [];

  const idGColListInput = ui.input.multiChoice('Seq id columns', {
    value: idDefaultValue,
    items: idGColList.map((gcol: DG.GridColumn) => gcol.name)
  });

  const seqGColList: DG.GridColumn[] = wu.count(0).take(grid.columns.length)/* range rom 0 to grid.columns.length */
    .map((colI: number) => grid.columns.byIndex(colI)!)
    .filter((gc: DG.GridColumn) => {
      const col: DG.Column | null = gc.column;
      if (col && col.semType === DG.SEMTYPE.MACROMOLECULE) {
        const sh = _package.seqHelper.getSeqHandler(col);
        return sh.isFasta();
      }
      return false;
    }).toArray();

  const seqDefaultValue = seqGColList.length > 0 ? seqGColList[0].name : [];
  const seqColInput = ui.input.choice('Seq column', {
    value: seqDefaultValue,
    items: seqGColList.map((gCol: DG.GridColumn) => gCol.name)
  });

  const lineWidthInput = ui.input.int('FASTA line width', {value: FASTA_LINE_WIDTH});

  ui.dialog({title: 'Save as FASTA'})
    .add(ui.inputs([
      idGColListInput,
      seqColInput,
      lineWidthInput,
    ]))
    .onOK(() => {
      const valueIdColList: DG.Column[] = idGColListInput.value ?
        idGColListInput.value.map((colName: string) => grid.columns.byName(colName)!.column!) : [];
      const valueSeqCol: DG.Column | null = seqColInput.value ?
        grid.columns.byName(seqColInput.value as string)!.column : null;
      const valueLineWidth = lineWidthInput.value ?? FASTA_LINE_WIDTH;

      if (!valueSeqCol)
        grok.shell.warning(`Seq column is mandatory to save as FASTA.`);

      const seqHandler = _package.seqHelper.getSeqHandler(valueSeqCol!);
      const resFastaTxt: string = saveAsFastaDo(valueIdColList, seqHandler, valueLineWidth);

      const aEl: HTMLAnchorElement = document.createElement('a');
      aEl.setAttribute('href', `data:text/plain;charset=utf-8,${encodeURIComponent(resFastaTxt)}`);
      aEl.setAttribute('download', `${grid.dataFrame.name}.fasta`);
      aEl.click();
    })
    .show();
}

/**
 * Builds FASTA content from id columns list and seq column
 * @param {DG.Column[]} idColList - list of columns with identifiers
 * @param {DG.Column} seqCol - column with sequence
 * @param {number} lineWidth - FASTA line width
 * @param {string} lineSeparator - FASTA line separator
 * @return {string} FASTA content
 */
export function saveAsFastaDo(
  idColList: DG.Column[], seqHandler: ISeqHandler, lineWidth: number = FASTA_LINE_WIDTH, lineSeparator: string = '\n',
): string {
  const fastaLines: string[] = [];

  for (let rowIdx: number = 0; rowIdx < seqHandler.length; rowIdx++) {
    // multiple identifiers separated by vertical bars
    // https://en.wikipedia.org/wiki/FASTA_format

    const seqId: string = idColList.map((col) => col.get(rowIdx).toString()).join('|');
    const srcSS = seqHandler.getSplitted(rowIdx);
    const seqLineList: string[] = wrapSequence(srcSS, lineWidth);

    fastaLines.push(`>${seqId}${lineSeparator}`);
    for (const line of seqLineList)
      fastaLines.push(`${line}${lineSeparator}`);
  }

  //return fastaLines.join(lineSeparator);
  return ''.concat(...fastaLines);
}

/* split sequence for monomers to prevent wrapping monomer partially */
export function wrapSequence(srcSS: ISeqSplitted, lineWidth: number = FASTA_LINE_WIDTH): string[] {
  let seqPos: number = 0;
  const seqLength: number = srcSS.length;

  const seqLineList: string[] = [];
  while (seqPos < seqLength) {
    /* join sliced monomer into line */
    const seqLine = wu.count(seqPos).take(Math.min(srcSS.length - seqPos, lineWidth)).map((p) => srcSS.getOriginal(p)).toArray();
    const seqLineTxt: string = seqLine.map((om) => om.length > 1 ? `[${om}]` : om)
      .reduce((a, b) => a + b, '');
    seqLineList.push(seqLineTxt);
    seqPos += seqLine.length;
  }

  return seqLineList;
}
