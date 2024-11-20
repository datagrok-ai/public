import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, expectObject}
  from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {doPolyToolConvert} from '../polytool/conversion/pt-conversion';
import {getOverriddenLibrary} from '../polytool/conversion/pt-synthetic';
import {getRules} from '../polytool/conversion/pt-rules';


import {_package} from '../package-test';

category('PolyTool: Convert', () => {
  let helmHelper: IHelmHelper;
  let seqHelper: ISeqHelper;
  let rdKitModule: RDModule;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; //backup

  before(async () => {
    helmHelper = await getHelmHelper();
    seqHelper = await getSeqHelper();
    rdKitModule = await getRdKitModule();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const tests = {
    'cyclized-C(1)-2-1': {
      src: {seq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)-meI'},
      tgt: {
        helm: 'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$V2.0',
        mol: {atomCount: 95, bondCount: 100, inchiKey: 'LMJUFVBPWWJJPN-AJJYTACESA-N'},
      }
    },
    'cyclized-C(1)-0-1': {
      src: {seq: 'C(1)-T-G-H-F-Y-P-C(1)-meI'},
      tgt: {
        helm: 'PEPTIDE1{C.T.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$V2.0',
        mol: {atomCount: 73, bondCount: 77, inchiKey: 'KLFRBMUPPMMGJM-HXTBFBBASA-N'},
      }
    },
    'cyclized-C(1)-2-0': {
      src: {seq: 'R-F-C(1)-T-G-H-F-Y-P-C(1)'},
      tgt: {
        helm: 'PEPTIDE1{R.F.C.T.G.H.F.Y.P.C}$PEPTIDE1,PEPTIDE1,3:R3-10:R3$$$V2.0',
        mol: {atomCount: 86, bondCount: 91, inchiKey: 'WIHSRTQGMICACU-DDDKLKPZSA-N'},
      }
    },
    'cyclized-C(1)-0-0': {
      src: {seq: 'C(1)-T-G-H-F-Y-P-C(1)'},
      tgt: {
        helm: 'PEPTIDE1{C.T.G.H.F.Y.P.C}$PEPTIDE1,PEPTIDE1,1:R3-8:R3$$$V2.0',
        mol: {atomCount: 64, bondCount: 68, inchiKey: 'LOSMDBLEXLWPLB-OFZKBENXSA-N'},
      }
    },
    'cyclized-D(2)-NH2(2)-3-0': {
      src: {seq: 'R-F-D(2)-T-G-H-F-Y-P-NH2(2)'},
      tgt: {
        helm: 'PEPTIDE1{R.F.D.T.G.H.F.Y.P.[NH2]}$PEPTIDE1,PEPTIDE1,10:R2-3:R3$$$V2.0',
        mol: {atomCount: 81, bondCount: 86, inchiKey: 'CBMGNYKOZWNVNK-AHGCAHLCSA-N'},
      }
    },
    'cyclized-D(2)-NH2(2)-0-0': {
      src: {seq: 'D(2)-T-G-H-F-Y-P-NH2(2)'},
      tgt: {
        helm: 'PEPTIDE1{D.T.G.H.F.Y.P.[NH2]}$PEPTIDE1,PEPTIDE1,8:R2-1:R3$$$V2.0',
        mol: {atomCount: 59, bondCount: 63, inchiKey: 'HGRHAUQBJXFERJ-MUFWPYSASA-N'},
      }
    },
    // 'cyclized-azG(4)-aG(4)-2-1': {
    //   src: {seq: 'R-F-azG(4)-T-G-H-F-Y-P-aG(4)-meI'},
    //   tgt: {
    //     helm: 'PEPTIDE1{R.F.[azG_GGaz].T.G.H.F.Y.P.[aG_GGaz].[meI]}|PEPTIDE2{[GGaz]}$PEPTIDE1,PEPTIDE2,3:R3-1:R1|PEPTIDE1,PEPTIDE2,10:R3-1:R2$$$V2.0',
    //     mol: {atomCount: 97, bondCount: 103, inchiKey: 'WJSYGVBGPCCSJF-PERUNASMSA-N'},
    //   }
    // },
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`toHelm-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const res = doPolyToolConvert([testData.src.seq], rules, helmHelper);
      expect(res[0], testData.tgt.helm);
    });
  }

  for (const [testName, testData] of Object.entries(tests)) {
    test(`toAtomicLevel-${testName}`, async () => {
      const rules = await getRules(['rules_example.json']);
      const helmList = doPolyToolConvert([testData.src.seq], rules, helmHelper);

      const lib = await getOverriddenLibrary(rules);

      const helmCol = DG.Column.fromStrings('helm', helmList);
      helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
      helmCol.meta.units = NOTATION.HELM;
      const talRes = await seqHelper.helmToAtomicLevel(helmCol, false, false, lib);
      if (talRes.warnings && talRes.warnings.length > 0)
        throw new Error(talRes.warnings[0]);
      expect(talRes.molCol != null, true, '.molCol is not null');
      const k = 11;
      const molfile = talRes.molCol!.get(0)!;
      const mol = rdKitModule.get_mol(molfile);
      try {
        const resMol = {
          atomCount: mol.get_num_atoms(), bondCount: mol.get_num_bonds(),
          inchiKey: rdKitModule.get_inchikey_for_inchi(mol.get_inchi())
        };
        expectObject(resMol, testData.tgt.mol);
      } finally {
        mol.delete();
      }
    });
  }

  test('ui-col-wo-table', async () => {
    const packageName = _package.name;
    const seqList = Object.values(tests).map((td) => td.src.seq);
    const helmList = Object.values(tests).map((td) => td.tgt.helm);
    const seqCol = DG.Column.fromStrings('seq', seqList);
    const helmCol: DG.Column = await grok.functions.call(`${packageName}:polyToolConvert2`, {
      seqCol: seqCol,
      generateHelm: true,
      chiralityEngine: true,
      rules: ['rules_example.json']
    });
    expectArray(helmCol.toList(), helmList);
    expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helmCol.meta.units, NOTATION.HELM);
    expect(helmCol.dataFrame, null);
  });

  test('ui-col', async () => {
    const packageName = _package.name;
    const seqList = Object.values(tests).map((td) => td.src.seq);
    const helmList = Object.values(tests).map((td) => td.tgt.helm);
    const seqCol = DG.Column.fromStrings('seq', seqList);
    const df = DG.DataFrame.fromColumns([seqCol]);
    const tv = grok.shell.addTableView(df);

    const helmCol: DG.Column = await grok.functions.call(`${packageName}:polyToolConvert2`, {
      seqCol: seqCol,
      generateHelm: true,
      chiralityEngine: true,
      rules: ['rules_example.json']
    });
    expectArray(helmCol.toList(), helmList);
    expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(helmCol.meta.units, NOTATION.HELM);
    expect(helmCol.dataFrame.id, df.id);

    await testEvent(tv.grid.onAfterDrawContent, () => {}, async () => {
      tv.grid.invalidate();
    }, 15000);
    expect(helmCol.getTag(DG.TAGS.CELL_RENDERER), 'helm');
  });
});
