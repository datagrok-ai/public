import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {TAGS as bioTAGS, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GAP_SYMBOL} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';

import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {buildCompositionTable} from '@datagrok-libraries/bio/src/utils/composition-table';
import '../../css/composition-analysis.css';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';

export function getCompositionAnalysisWidget(
  val: DG.SemanticValue, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper
): DG.Widget {
  const host = ui.div();
  host.classList.add('macromolecule-cell-comp-analysis-host');
  const alphabet = val.cell.column.tags[bioTAGS.alphabet];
  const biotype = alphabet === ALPHABET.DNA || alphabet === ALPHABET.RNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;

  const counts: { [m: string]: number } = {};
  const sh = seqHelper.getSeqHandler(val.cell.column as DG.Column<string>);
  const rowIdx = val.cell.rowIndex;
  const seqSS = sh.getSplitted(rowIdx);
  // in case of HELM, there might be multiple biotypes in one sequence
  const bioTypes: {[symbol: string]: HelmType} = {};
  wu.count(0).take(seqSS.length).filter((posIdx) => !seqSS.isGap(posIdx)).forEach((posIdx) => {
    let cm = seqSS.getCanonical(posIdx);
    if (biotype === HelmTypes.NUCLEOTIDE && sh.isHelm() && cm[1] === '(' && cm[cm.length - 2] === ')')
      cm = cm.substring(2, cm.length - 2);
    const count = counts[cm] || 0;
    counts[cm] = count + 1;
    if (!bioTypes[cm] && seqSS.graphInfo?.polymerTypes) {
      const polymerType = seqSS.graphInfo.polymerTypes[posIdx];
      bioTypes[cm] = polymerTypeToHelmType(polymerType);
    }
  });


  const table = buildCompositionTable(counts, biotype, monomerLib, Object.keys(bioTypes).length ? bioTypes : undefined);
  Array.from(table.rows).forEach((row) => {
    const barCol = (row.getElementsByClassName('macromolecule-cell-comp-analysis-bar')[0] as HTMLDivElement)
      .style.backgroundColor;
    row.cells[0].style.color = barCol;
  });

  host.appendChild(table);
  return new DG.Widget(host);
}

