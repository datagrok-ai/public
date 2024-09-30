import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {TAGS as bioTAGS, ALPHABET, getPaletteByType} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import '../../css/composition-analysis.css';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {GAP_SYMBOL} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';


export function getCompositionAnalysisWidget(val: DG.SemanticValue): DG.Widget {
  const host = ui.div();
  host.classList.add('macromolecule-cell-comp-analysis-host');
  const alphabet = val.cell.column.tags[bioTAGS.alphabet];
  let palette: SeqPalette = UnknownSeqPalettes.Color;
  switch (alphabet) {
  case ALPHABET.DNA:
  case ALPHABET.RNA:
    palette = getPaletteByType(ALPHABET.DNA);
    break;
  case ALPHABET.PT:
    palette = getPaletteByType(ALPHABET.PT);
    break;
  default:
    break;
  }

  const counts: { [m: string]: number } = {};
  const sh = SeqHandler.forColumn(val.cell.column as DG.Column<string>);
  const rowIdx = val.cell.rowIndex;
  const seqSS = sh.getSplitted(rowIdx);
  wu.count(0).take(seqSS.length).filter((posIdx) => !seqSS.isGap(posIdx)).forEach((posIdx) => {
    const cm = seqSS.getCanonical(posIdx);
    const count = counts[cm] || 0;
    counts[cm] = count + 1;
  });
  const table = buildCompositionTable(palette, counts);
  Array.from(table.rows).forEach((row) => {
    const barCol = (row.getElementsByClassName('macromolecule-cell-comp-analysis-bar')[0] as HTMLDivElement)
      .style.backgroundColor;
    row.cells[0].style.color = barCol;
  });

  host.appendChild(table);
  return new DG.Widget(host);
}

export function buildCompositionTable(palette: SeqPalette, counts: { [m: string]: number }): HTMLTableElement {
  let sumValue: number = 0;
  let maxValue: number | null = null;
  for (const value of Object.values(counts)) {
    sumValue = sumValue + value;
    maxValue = maxValue === null ? value : Math.max(maxValue, value);
  }
  const maxRatio = maxValue! / sumValue;
  const elMap: { [m: string]: HTMLElement } = Object.assign({}, ...Array.from(Object.entries(counts))
    .sort((a, b) => b[1] - a[1])
    .map(([cm, value]) => {
      const ratio = value / sumValue;
      const color = palette.get(cm);
      const barDiv = ui.div('', {classes: 'macromolecule-cell-comp-analysis-bar'});
      barDiv.style.width = `${50 * ratio / maxRatio}px`;
      barDiv.style.backgroundColor = color;
      if (GAP_SYMBOL === cm) {
        barDiv.style.borderWidth = '1px';
        barDiv.style.borderStyle = 'solid';
        barDiv.style.borderColor = DG.Color.toHtml(DG.Color.lightGray);
      }
      const displayMonomer: string = GAP_SYMBOL === cm ? '-' : cm;
      const valueDiv = ui.div(`${(100 * ratio).toFixed(2)}%`);
      const el = ui.div([barDiv, valueDiv], {classes: 'macromolecule-cell-comp-analysis-value'});
      return ({[displayMonomer]: el});
    }));

  const table = ui.tableFromMap(elMap);
  Array.from(table.rows).forEach((row) => {
    const barCol = (row.getElementsByClassName('macromolecule-cell-comp-analysis-bar')[0] as HTMLDivElement)
      .style.backgroundColor;
    row.cells[0].style.color = barCol;
  });
  return table;
}
