import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TAGS as bioTAGS, ALPHABET} from './macromolecule';
import {GAP_SYMBOL} from './macromolecule/consts';
import {IMonomerLibBase} from './../types/monomer-library';
import {HelmType} from './../helm/types';
import {helmTypeToPolymerType} from '../monomer-works/monomer-works';

export function buildCompositionTable(
  counts: { [m: string]: number }, biotype: HelmType, monomerLib: IMonomerLibBase,
  monomerBioTypes?: { [m: string]: HelmType }
): HTMLTableElement {
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
      let color: string;
      const actBioType = monomerBioTypes ? (monomerBioTypes[cm] || biotype) : biotype;
      try {
        const colors = monomerLib.getMonomerColors(actBioType, cm);
        color = colors?.backgroundcolor || '#CCCCCC';
      } catch (error) {
        console.warn(`Failed to get colors for monomer ${cm}:`, error);
        color = '#CCCCCC';
      }

      const barDiv = ui.div('', {classes: 'macromolecule-cell-comp-analysis-bar'});
      barDiv.style.width = `${50 * ratio / maxRatio}px`;
      barDiv.style.backgroundColor = color;
      if (GAP_SYMBOL === cm) {
        barDiv.style.borderWidth = '1px';
        barDiv.style.borderStyle = 'solid';
        barDiv.style.borderColor = DG.Color.toHtml(DG.Color.lightGray);
      }
      let monomerDisplayName = cm;
      if (cm !== GAP_SYMBOL) {
        const m = monomerLib.getMonomer(helmTypeToPolymerType(actBioType), cm);
        if (m && m.symbol !== cm) // for explicit smiles based monomers
          monomerDisplayName = m.symbol;
      } else
        monomerDisplayName = '-';
      const valueDiv = ui.div(`${(100 * ratio).toFixed(2)}%`);
      const el = ui.div([barDiv, valueDiv], {classes: 'macromolecule-cell-comp-analysis-value'});
      return ({[monomerDisplayName]: el});
    }));

  const table = ui.tableFromMap(elMap);
  Array.from(table.rows).forEach((row) => {
    const barCol = (row.getElementsByClassName('macromolecule-cell-comp-analysis-bar')[0] as HTMLDivElement)
      .style.backgroundColor;
    row.cells[0].style.color = barCol;
  });
  return table;
}
