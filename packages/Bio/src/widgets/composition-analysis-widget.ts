import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {TAGS as bioTAGS, ALPHABET, getPaletteByType} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import '../../css/composition-analysis.css';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';


export function getCompositionAnalysisWidget(val: DG.SemanticValue) {
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

  const counts = new Map<string, number>();
  let max = 0;
  const uh = UnitsHandler.getOrCreate(val.cell.column);
  const splitter = uh.getSplitter();
  const parts = splitter(val.value);
  const len = parts.length;
  wu(parts).filter((p) => !!p && p !== '').forEach((c: string) => {
    const count = counts.get(c) || 0;
    counts.set(c, count + 1);
    max = Math.max(max, count + 1);
  });
  max /= len;// percentage
  // calculate frequencies
  const compositionMap: { [key: string]: HTMLElement } = {};
  const valueArray = Array.from(counts.entries());
  valueArray.sort((a, b) => b[1] - a[1]);
  valueArray.forEach(([key, value]) => {
    const ratio = value / len;
    const color = palette.get(key);
    const barDiv = ui.div('', {classes: 'macromolecule-cell-comp-analysis-bar'});
    barDiv.style.width = `${ratio / max * 50}px`;
    const valueDiv = ui.div((ratio * 100).toFixed(2) + '%');
    barDiv.style.backgroundColor = color;
    compositionMap[key] = ui.div([barDiv, valueDiv], {classes: 'macromolecule-cell-comp-analysis-value'});
  });

  const table = ui.tableFromMap(compositionMap);
  Array.from(table.rows).forEach((row) => {
    const barCol = (row.getElementsByClassName('macromolecule-cell-comp-analysis-bar')[0] as HTMLDivElement)
      .style.backgroundColor;
    row.cells[0].style.color = barCol;
  });

  host.appendChild(table);
  return new DG.Widget(host);
}
