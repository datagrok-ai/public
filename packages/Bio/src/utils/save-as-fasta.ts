import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import wu from 'wu';
import {splitterAsFasta, SplitterFunc, UnitsHandler} from '@datagrok-libraries/bio';

const FASTA_LINE_WIDTH = 60;

/** Shows dialog to select id columns list and seq column, builds and downloads FASTA content */
export function saveAsFastaUI() {
  // Use grid for column order adjusted by user
  const grid: DG.Grid = grok.shell.tv.grid;

  const idGColList: DG.GridColumn[] = wu.count(0).take(grid.columns.length)
    .map((colI: number) => grid.columns.byIndex(colI)!)
    .filter((gcol: DG.GridColumn) => gcol.column ? gcol.column.semType !== DG.SEMTYPE.MACROMOLECULE : false).toArray();
  const defaultIdGCol: DG.GridColumn | undefined = idGColList
    .find((gcol: DG.GridColumn) => gcol.name.toLowerCase().indexOf('id') !== -1);
  const idDefaultValue = defaultIdGCol ? [defaultIdGCol.name] : [];

  const idGColListInput = ui.multiChoiceInput('Seq id columns', idDefaultValue,
    idGColList.map((gcol: DG.GridColumn) => gcol.name));

  const seqGColList: DG.GridColumn[] = wu.count(0).take(grid.columns.length)/* range rom 0 to grid.columns.length */
    .map((colI: number) => grid.columns.byIndex(colI)!)
    .filter((gc: DG.GridColumn) => {
      const col: DG.Column | null = gc.column;
      if (col && col.semType === DG.SEMTYPE.MACROMOLECULE) {
        const uh = new UnitsHandler(col);
        return uh.isFasta();
      }
      return false;
    }).toArray();

  const seqDefaultValue = seqGColList.length > 0 ? seqGColList[0].name : [];
  const seqColInput = ui.choiceInput('Seq column', seqDefaultValue,
    seqGColList.map((gCol: DG.GridColumn) => gCol.name));

  const lineWidthInput = ui.intInput('FASTA line width', FASTA_LINE_WIDTH);

  ui.dialog({title: 'Save as FASTA'})
    .add(ui.inputs([
      idGColListInput,
      seqColInput,
      lineWidthInput
    ]))
    .onOK(() => {
      const valueIdColList: DG.Column[] = idGColListInput.value ?
        idGColListInput.value.map((colName: string) => grid.columns.byName(colName)!.column!) : [];
      const valueSeqCol: DG.Column | null = seqColInput.value ?
        grid.columns.byName(seqColInput.value as string)!.column : null;
      const valueLineWidth = lineWidthInput.value ?? FASTA_LINE_WIDTH;

      if (!valueSeqCol)
        grok.shell.warning(`Seq column is mandatory to save as FASTA.`);

      const resFastaTxt: string = saveAsFastaDo(valueIdColList, valueSeqCol!, valueLineWidth);

      const aEl: HTMLAnchorElement = document.createElement('a');
      aEl.setAttribute('href', `data:text/plain;charset=utf-8,${encodeURIComponent(resFastaTxt)}`);
      aEl.setAttribute('download', `${grid.dataFrame.name}.fasta`);
      aEl.click();
    })
    .show();
}

/** */
export function saveAsFastaDo(
  idColList: DG.Column[], seqCol: DG.Column, lineWidth: number = FASTA_LINE_WIDTH, lineSeparator: string = '\n'
): string {
  const splitter: SplitterFunc = splitterAsFasta;

  const fastaLines: string[] = [];

  for (let rowI: number = 0; rowI < seqCol.length; rowI++) {
    // multiple identifiers separated by vertical bars
    // https://en.wikipedia.org/wiki/FASTA_format

    const seqId: string = idColList.map((col) => col.get(rowI).toString()).join('|');
    const seq: string = seqCol.get(rowI);
    const seqLineList: string[] = wrapSequence(seq, splitter, lineWidth);

    fastaLines.push(`>${seqId}${lineSeparator}`);
    for (const line of seqLineList)
      fastaLines.push(`${line}${lineSeparator}`);
  }

  //return fastaLines.join(lineSeparator);
  return ''.concat(...fastaLines);
}

/* split sequence for monomers to prevent wrapping monomer partially */
export function wrapSequence(seq: string, splitter: SplitterFunc, lineWidth: number = FASTA_LINE_WIDTH): string[] {
  const seqMonomerList = splitter(seq);
  let seqPos: number = 0;
  const seqLength: number = seqMonomerList.length;

  const seqLineList: string[] = [];
  while (seqPos < seqLength) {
    /* join sliced monomer into line */
    const seqLine: string[] = seqMonomerList.slice(seqPos, seqPos + lineWidth);
    const seqLineTxt: string = seqLine.map((m) => m.length > 1 ? `[${m}]` : m).join('');
    seqLineList.push(seqLineTxt);
    seqPos += seqLine.length;
  }

  return seqLineList;
}
