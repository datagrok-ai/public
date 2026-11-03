import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getMonomerLibHelper, IMonomerLib, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {assembleRowWithMap, assembleSeqColumn} from '@datagrok-libraries/bio/src/monomer-works/assembler';

/** The monomer libraries Bio registers with Chem's macromolecule assembler, which keeps only the last few: a key it
 * has dropped (every PolyTool conversion and rule preview makes an override library of its own) is registered again,
 * not refused with "no monomer library ...: setLibrary first". */
category('toAtomicLevelCrux: libraries', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let seqHelper: ISeqHelper;
  let monomerLib: IMonomerLib;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests(); // HELMCoreLibrary.json + polytool-lib.json
    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  test('a library the assembler dropped is registered again', async () => {
    const col = DG.Column.fromList(DG.TYPE.STRING, 'seq', ['PEPTIDE1{A.G.C}$$$$V2.0']) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = NOTATION.HELM;
    // an override library of its own per conversion, as PolyTool makes them: more than the assembler keeps (8)
    const first = monomerLib.override({}, 'toAtomicLevelCrux: libraries, first');
    const firstRow = await assembleRowWithMap(col, 0, first, seqHelper, 'linear');
    expect(firstRow.error, null, `the first library: ${firstRow.error}`);
    for (let i = 0; i < 10; i++) {
      const row = await assembleRowWithMap(col, 0, monomerLib.override({}, `toAtomicLevelCrux: libraries, ${i}`),
        seqHelper, 'linear');
      expect(row.error, null, `library ${i}: ${row.error}`);
    }
    // the first library's column again: its hover maps and a new conversion
    const again = await assembleRowWithMap(col, 0, first, seqHelper, 'linear');
    expect(again.error, null, `the first library again: ${again.error}`);
    expect(again.molfile, firstRow.molfile, 'the same molecule');
    expect(again.monomers.size, 3, 'the monomer map');
    const res = await assembleSeqColumn(col, first, seqHelper, 'linear', false);
    expect(res.warnings.length, 0, res.warnings.join('; '));
    expect(res.molCol?.get(0), firstRow.molfile, 'the column\'s molecule');
  }, {timeout: 60000});
});
