import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import * as ngl from 'NGL';

import {category, expect/*, expect*/, expectObject, test} from '@datagrok-libraries/utils/src/test';
import {Molecule3DUnitsHandler} from '@datagrok-libraries/bio/src/molecule-3d';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {IPdbAtomBase, IPdbqtAtomBase} from '@datagrok-libraries/bio/src/pdb/format/types';
import {AtomBase, AtomCoordsBase, LineBase} from '@datagrok-libraries/bio/src/pdb/format/types-base';
import {PdbAtomCoords} from '@datagrok-libraries/bio/src/pdb/format/types-pdb';

import {IMPORT} from '../consts-import';
import {importPdbqtUI} from '../utils/import-pdbqt';
import {PdbqtAtomCoords} from '@datagrok-libraries/bio/src/pdb/format/types-pdbqt';
import {Pdbqt} from '../utils/pdbqt-parser';

import {_package} from '../package-test';

category('pdbqt', () => {
  const pdbqtAtomTestData: {
    [name: string]: { pdbqtStr: string, pdbqtAtom: IPdbqtAtomBase, pdbAtom: IPdbAtomBase, pdbStr: string }
  } = {
    simple: {
      /* atom line from ligand_out.pdbqt was fixed according to 23-26 positions for resNum aligned right */
      pdbqtStr: /**/ 'ATOM      1  C   LIG     1      16.882  -0.414  35.505  0.00  0.00    -0.069 C ',
      pdbStr: /*  */ 'ATOM      1  C   LIG     1      16.882  -0.414  35.505  0.00  0.00           C  ',
      // @formatter:off
      pdbqtAtom: new PdbqtAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        1, 'C', '', '', 'LIG', '', 1, ''),
      16.882, -0.414, 35.505, 0, 0),
      -0.069, 'C'),
      pdbAtom: new PdbAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        1, 'C', '', '', 'LIG', '', 1, ''),
      16.882, -0.414, 35.505, 0, 0),
      '', 'C', ''),
      // @formatter:on
    },
    alpha: {
      /* atom line from ligand_out.pdbqt was fixed according to 23-26 positions for resNum aligned right */
      pdbqtStr: /**/ 'ATOM     20  O   LIG     1      15.199  -3.634  34.355  0.00  0.00    -0.390 OA',
      pdbStr: /*  */ 'ATOM     20  O   LIG     1      15.199  -3.634  34.355  0.00  0.00           O  ', // from NGL
      // @formatter:off
      pdbqtAtom: new PdbqtAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        20, 'O', '', '', 'LIG', '', 1, ''),
      15.199, -3.634, 34.355, 0, 0),
      -0.390, 'OA'),
      pdbAtom: new PdbAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        20, 'O', '', '', 'LIG', '', 1, ''),
      15.199, -3.634, 34.355, 0, 0),
      '', 'O', ''),
      // @formatter:on
    },
    hydrogen: {
      /* atom line from 3SWZ.pdbqt was fixed for partial_charge to have explicit sign*/
      pdbqtStr: /**/ 'ATOM      7  HN2 SER A  30      69.016  10.250  16.854  0.00  0.00    +0.139 HD',
      pdbStr: /*  */ 'ATOM      7  HN2 SER A  30      69.016  10.250  16.854  0.00  0.00           H  ', // from NGL
      // @formatter:off
      pdbqtAtom: new PdbqtAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        7, 'H', 'N2', '', 'SER', 'A', 30, ''),
      69.016, 10.250, 16.854, 0, 0),
      0.139, 'HD'),
      pdbAtom: new PdbAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        7, 'H', 'N2', '', 'SER', 'A', 30, ''),
      69.016, 10.250, 16.854, 0, 0),
      '', 'H', ''),
      // @formatter:on
    },
  };

  for (const [testName, testData] of Object.entries(pdbqtAtomTestData)) {
    test(`${testName}-pdbqtAtomFromStr`, async () => {
      const atomRes = PdbqtAtomCoords.fromStr(testData.pdbqtStr);
      expectObject(atomRes, testData.pdbqtAtom);
    });

    test(`${testName}-pdbqtAtomToStr`, async () => {
      const pdbqtStrRes = testData.pdbqtAtom.toStr();
      expect(pdbqtStrRes, testData.pdbqtStr);
    });

    test(`${testName}-pdbqtAtomToPdb`, async () => {
      const pdbAtomRes: IPdbAtomBase = testData.pdbqtAtom.toPdb();
      expectObject(pdbAtomRes, testData.pdbAtom);
    });

    test(`${testName}-pdbAtomToStr`, async () => {
      const pdbStrRes = testData.pdbAtom.toStr();
      expect(pdbStrRes, testData.pdbStr);
    });

    // // Not implemented
    // test(`pdbAtomFromStr-${testName}`, async () => {
    //   const atomRes = PdbAtom.fromStr(testData.pdbStr);
    //   expectObject(atomRes, testData.pdbAtom);
    // });
  }

  test('parse', async () => {
    await _testPdbqtParse();
  });

  test('parse-autodock-gpu', async () => {
    await _testPdbqtParseAutodockGpu();
  });

  test('import-models', async () => {
    await _testPdbqtImportModels();
  });

  test('import-target-only', async () => {
    await _testPdbqtImportTargetOnly();
  });

  // -- NGL PdbqtParser --

  test('ngl-pdbqt-parser-structure', async () => {
    await _testNglPdbqtParserStructure();
  });
  test('ngl-pdbqt-parser-ligands', async () => {
    await _testNglPdbqtParserLigands();
  });
});

async function _testPdbqtParse(): Promise<void> {
  try {
    const cnt: string = await _package.files.readAsText('docking/ligand_out.pdbqt');
    const data: Pdbqt = Pdbqt.parse(cnt);

    expect(data['currentModel'], null);

    expect(data.models.length, 10);
    for (const model of data.models) {
      expect(model.torsdof, 3);
      expect(model.atoms.length, 19);
      expect(model.children.length, 2);
      expect(model.children[0].atoms.length, 2);
      expect(model.children[0].children.length, 0);

      expect(model.children[1].atoms.length, 1);
      expect(model.children[1].children.length, 1);
      expect(model.children[1].children[0].atoms.length, 8);
      expect(model.children[1].children[0].children.length, 0);

      const flatAtomList: PdbqtAtomCoords[] = [];
      model.flattenAtoms(flatAtomList);
      expect(flatAtomList.length, 30);
    }
  } catch (err) {
    const [errMsg, errStack] = errInfo(err);
    _package.logger.error(errMsg, undefined, errStack);
    throw err;
  }
}

async function _testPdbqtParseAutodockGpu(): Promise<void> {
  const cnt = await _package.files.readAsText('samples/1bdq.autodock-gpu.pdbqt');
  const data: Pdbqt = Pdbqt.parse(cnt);

  expect(data['currentModel'], null); // last model closed
  expect(data.models.length, 2);
  // for (const model of data.models) {
  //   expect(model.torsdof, 14);
  //   expect(model.atoms.length, 19);
  //   expect(model.children.length, 2);
  // }
}

async function _testPdbqtImportModels(): Promise<void> {
  try {
    const cnt: string = await _package.files.readAsText('docking/ligand_out.pdbqt');
    const df: DG.DataFrame = (await importPdbqtUI(cnt, true))[0];

    const molCol = df.getCol(IMPORT.pdb.molColName);
    const uh = IMPORT.pdb.unitsHandlerClass.getOrCreate(molCol);
    expect(uh.units, IMPORT.pdb.units);
  } catch (err) {
    const [errMsg, errStack] = errInfo(err);
    _package.logger.error(errMsg, undefined, errStack);
    throw err;
  }
}

async function _testPdbqtImportTargetOnly(): Promise<void> {
  const cnt: string = await _package.files.readAsText('docking/3SWZ.pdbqt');
  const dfList = await importPdbqtUI(cnt, true);
  expect(dfList.length, 0);

  const view = grok.shell.v;
  expect(grok.shell.v.name, 'Mol*');
}

async function _testNglPdbqtParserStructure(): Promise<void> {
  const cntStr = await _package.files.readAsText('docking/3SWZ.pdbqt');
  const cntBlob = new Blob([cntStr]);

  // from NGL.autoLoad()
  const val = await ngl.autoLoad(cntBlob, {ext: 'pdbqt'});
  const pdbStr = (new ngl.PdbWriter(val)).getString();

  expect(val.atomCount, 4651);
  expect(val.bondCount, 4913);
}

async function _testNglPdbqtParserLigands(): Promise<void> {
  const cntStr = await _package.files.readAsText('docking/ligand_out.pdbqt');
  const cntBlob = new Blob([cntStr]);

  // from NGL.autoLoad()
  const val: ngl.Structure = await ngl.autoLoad(cntBlob, {ext: 'pdbqt'});

  expect(val.modelStore.count, 10);
  const pose0Val = val.getView(new ngl.Selection('/0'));
  expect(pose0Val.atomCount, 30);
  expect(pose0Val.bondCount, 34);
}
