import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';

import {getCompositionAnalysisWidget} from '../widgets/composition-analysis-widget';

category('compositionAnalysis', () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
  });

  test('moreThan20Monomers', async () => {
    const df = DG.DataFrame.fromCsv(`seq\nACDEFGHIKLMNPQRSTVWYX`);
    const seqCol = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.UN);

    const widget = getCompositionAnalysisWidget(DG.SemanticValue.fromTableCell(df.cell(0, 'seq')),
      monomerLibHelper.getMonomerLib(), seqHelper);
    const table = widget.root.querySelector('table')!;
    expect(table.getElementsByClassName('macromolecule-cell-comp-analysis-bar').length, 20);
    expect(table.rows.length, 21);
  });
});
