import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {category, expect/*, expect*/, expectObject, test} from '@datagrok-libraries/utils/src/test';
import {Molecule3DUnitsHandler} from '@datagrok-libraries/bio/src/molecule-3d';

import {Pdbqt} from '../utils/pdbqt/pdbqt-parser';
import {errInfo} from '../utils/err-info';
import {importPdbqt} from '../package';
import {IMPORT} from '../consts-import';
import {IPdbAtomBase, IPdbqtAtomBase} from '../utils/pdbqt/types';
import {AtomBase, AtomCoordsBase, LineBase} from '../utils/pdbqt/types-base';
import {PdbqtAtomCoords} from '../utils/pdbqt/types-pdbqt';
import {PdbAtomCoords} from '../utils/pdbqt/types-pdb';

import {_package} from '../package-test';

category('pdbqt', () => {
  const pdbqtAtomTestData: {
    [name: string]: { pdbqtStr: string, pdbqtAtom: IPdbqtAtomBase, pdbAtom: IPdbAtomBase, pdbStr: string }
  } = {
    simple: {
      /* atom lines from ligand_out.pdbqt was fixed according to 23-26 positions for resNum aligned right */
      pdbqtStr: 'ATOM      1  C   LIG     1      16.882  -0.414  35.505  0.00  0.00    -0.069 C ',
      // @formatter:off
      pdbqtAtom: new PdbqtAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        1, 'C', '', '', 'LIG', '', 1, ''),
      16.882, -0.414, 35.505, 0, 0),
      -0.069, 'C'),
      // @formatter:on
      pdbStr: 'ATOM      1  C   LIG     1      16.882  -0.414  35.505  0.00  0.00           C  ',
      // @formatter:off
      pdbAtom: new PdbAtomCoords(new AtomCoordsBase(new AtomBase(new LineBase('ATOM'),
        1, 'C', '', '', 'LIG', '', 1, ''),
      16.882, -0.414, 35.505, 0, 0),
      '', 'C', ''),
      // @formatter:on
    },
  };

  for (const [testName, testData] of Object.entries(pdbqtAtomTestData)) {
    test(`pdbqtAtomFromStr-${testName}`, async () => {
      const atomRes = PdbqtAtomCoords.fromStr(testData.pdbqtStr);
      expectObject(atomRes, testData.pdbqtAtom);
    });

    test(`pdbqtAtomToStr-${testName}`, async () => {
      const pdbqtStrRes = testData.pdbqtAtom.toStr();
      expect(pdbqtStrRes, testData.pdbqtStr);
    });

    test(`pdbqtAtomToPdb-${testName}`, async () => {
      const pdbAtomRes: IPdbAtomBase = testData.pdbqtAtom.toPdb();
      expectObject(pdbAtomRes, testData.pdbAtom);
    });


    test(`pdbAtomToStr-${testName}`, async () => {
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
    await _testPdbqtParser();
  });

  // test('import-models', async () => {
  //   await _testPdbqtImportModels();
  // });

  test('import-target-only', async () => {
    await _testPdbqtImportTargetOnly();
  });
})
;

async function _testPdbqtParser(): Promise<void> {
  try {
    //const cnt: string = await _package.files.readAsText('docking/ligand_out.pdbqt');
    const dapiFilePath = `System:AppData/${_package.name}/docking/ligand_out.pdbqt`;
    const cnt = await grok.dapi.files.readAsText(dapiFilePath);
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

// async function _testPdbqtImportModels(): Promise<void> {
//   try {
//     const cnt: string = await _package.files.readAsText('docking/ligand_out.pdbqt');
//     const df: DG.DataFrame = (await importPdbqt(cnt))[0];
//
//     const molCol = df.getCol(IMPORT.pdb.molColName);
//     const uh = IMPORT.pdb.unitsHandlerClass.getOrCreate(molCol);
//     expect(uh.units, IMPORT.pdb.units);
//   } catch (err) {
//     const [errMsg, errStack] = errInfo(err);
//     _package.logger.error(errMsg, undefined, errStack);
//     throw err;
//   }
// }

async function _testPdbqtImportTargetOnly(): Promise<void> {
  const cnt: string = await _package.files.readAsText('docking/3SWZ.pdbqt');
  const dfList = await importPdbqt(cnt);
  expect(dfList.length, 0);

  const view = grok.shell.v;
  expect(grok.shell.v.name, 'Mol*');
}
