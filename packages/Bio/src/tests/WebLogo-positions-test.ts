import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test, testEvent} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {
  countForMonomerAtPosition,
  PositionInfo as PI,
  PositionMonomerInfo as PMI,
  WebLogoViewer,
} from '../viewers/web-logo-viewer';
import {GAP_SYMBOL} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';

const g: string = GAP_SYMBOL;

category('WebLogo-positions', () => {
  const csvDf1 = `seq
ATC-G-TTGC--
ATC-G-TTGC--
-TC-G-TTGC--
-TC-GCTTGC--
-TC-GCTTGC--`;

  test('allPositions', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf1);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ.MSA');

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo')) as WebLogoViewer;
    await testEvent(wlViewer.onLayoutCalculated, () => {}, () => {
      tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);
    }, 500, 'Layout calculate timeout');
    const positions: PI[] = wlViewer['positions'];

    const resAllDf1: PI[] = [
      new PI(0, '1', {'A': new PMI(2), [g]: new PMI(3)}),
      new PI(1, '2', {'T': new PMI(5)}),
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(3, '4', {[g]: new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {[g]: new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
      new PI(7, '8', {'T': new PMI(5)}),
      new PI(8, '9', {'G': new PMI(5)}),
      new PI(9, '10', {'C': new PMI(5)}),
      new PI(10, '11', {[g]: new PMI(5)}),
      new PI(11, '12', {[g]: new PMI(5)}),
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const m of positions[i].getMonomers())
        expect(positions[i].getFreq(m).rowCount, resAllDf1[i].getFreq(m).rowCount);
    }
    await wlViewer.awaitRendered();
  });

  test('positions with shrinkEmptyTail option true (filtered)', async () => {
    const csvDf2 = `seq
-TC-G-TTGC--
-TC-GCTTGC--
-T--C-GT-
-T--C-GT-
-T--C-GT-
-T--CCGT-`;
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf2);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ');

    df.filter.init((i) => {
      return i > 2;
    });
    df.filter.fireChanged();
    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'shrinkEmptyTail': true})) as unknown as WebLogoViewer;
    await testEvent(wlViewer.onLayoutCalculated, () => {}, () => {
      tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);
    }, 500);
    const positions: PI[] = wlViewer['positions'];

    const resAllDf1: PI[] = [
      new PI(0, '1', {[g]: new PMI(3)}),
      new PI(1, '2', {'T': new PMI(3)}),
      new PI(2, '3', {[g]: new PMI(3)}),
      new PI(3, '4', {[g]: new PMI(3)}),
      new PI(4, '5', {'C': new PMI(3)}),
      new PI(5, '6', {[g]: new PMI(2), 'C': new PMI(1)}),
      new PI(6, '7', {'G': new PMI(3)}),
      new PI(7, '8', {'T': new PMI(3)}),
      new PI(8, '9', {[g]: new PMI(3)}),
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const m of positions[i].getMonomers())
        expect(positions[i].getFreq(m).rowCount, resAllDf1[i].getFreq(m).rowCount);
    }
    await wlViewer.awaitRendered();
  });

  test('positions with skipEmptyPositions option', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf1);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ.MSA');

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'skipEmptyPositions': true})) as unknown as WebLogoViewer;
    await testEvent(wlViewer.onLayoutCalculated, () => {}, () => {
      tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);
    }, 500);
    const resPosList: PI[] = wlViewer['positions'];

    const tgtPosList: PI[] = [
      new PI(0, '1', {'A': new PMI(2), [g]: new PMI(3)}),
      new PI(1, '2', {'T': new PMI(5)}),
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {[g]: new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
      new PI(7, '8', {'T': new PMI(5)}),
      new PI(8, '9', {'G': new PMI(5)}),
      new PI(9, '10', {'C': new PMI(5)}),
    ];

    expect(resPosList.length, tgtPosList.length);
    for (let posI = 0; posI < resPosList.length; posI++) {
      const resPos = resPosList[posI];
      const tgtPos = tgtPosList[posI];
      expectPositionInfo(resPos, tgtPos);
    }
    await wlViewer.awaitRendered();
  });

  test('count sequences for monomer at position', async () => {
    const df: DG.DataFrame = buildDfWithSeqCol(csvDf1, NOTATION.FASTA, ALPHABET.DNA, 'SEQ.MSA');
    const seqCol: DG.Column = df.getCol('seq');

    const tv: DG.TableView = grok.shell.addTableView(df);

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo', {
      startPositionName: '3',
      endPositionName: '7',
      skipEmptyPositions: true,
    })) as unknown as WebLogoViewer;
    await testEvent(wlViewer.onLayoutCalculated, () => {}, () => {
      tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);
    }, 500);
    const resPosList: PI[] = wlViewer['positions'];
    const tgtPosList: PI[] = [
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {[g]: new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
    ];

    expect(resPosList.length, tgtPosList.length);
    for (let posI = 0; posI < resPosList.length; posI++) {
      const resPos = resPosList[posI];
      const tgtPos = tgtPosList[posI];
      expectPositionInfo(resPos, tgtPos);
    }

    const atPI1: PI = resPosList[1];
    const sh = SeqHandler.forColumn(seqCol);
    const countAt1 = countForMonomerAtPosition(df, sh, df.filter, 'G', atPI1);
    expect(countAt1, 5);
    await wlViewer.awaitRendered();
  });

  test('empty', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromColumns([(() => {
      const col = DG.Column.fromStrings('seq', []);
      col.semType = DG.SEMTYPE.MACROMOLECULE;
      col.meta.units = NOTATION.FASTA;
      col.setTag(bioTAGS.alphabet, ALPHABET.DNA);
      return col;
    })()]);

    const tv: DG.TableView = grok.shell.addTableView(df);

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo')) as WebLogoViewer;
    await testEvent(wlViewer.onLayoutCalculated, () => {}, () => {
      tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);
    }, 500);
    const resPosList: PI[] = wlViewer['positions'];
    await wlViewer.awaitRendered();
  });
});

function expectPositionInfo(actualPos: PI, expectedPos: PI): void {
  expect(actualPos.name, expectedPos.name);
  expectArray(actualPos.getMonomers(), expectedPos.getMonomers());
  for (const key of actualPos.getMonomers()) {
    //
    expect(actualPos.getFreq(key).rowCount, expectedPos.getFreq(key).rowCount);
  }
}

function buildDfWithSeqCol(csv: string, notation: NOTATION, alphabet: ALPHABET, aligned: string): DG.DataFrame {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);

  const seqCol: DG.Column = df.getCol('seq');
  seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  seqCol.meta.units = notation;
  seqCol.setTag(bioTAGS.alphabet, alphabet);
  seqCol.setTag(bioTAGS.aligned, aligned);

  return df;
}
