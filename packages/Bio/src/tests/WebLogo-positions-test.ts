import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, NOTATION, SplitterFunc, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  countForMonomerAtPosition,
  FilterSources,
  PositionInfo as PI,
  PositionMonomerInfo as PMI,
  WebLogoViewer
} from '../viewers/web-logo-viewer';

category('WebLogo-positions', () => {
  let tvList: DG.TableView[];
  let dfList: DG.DataFrame[];
  let currentView: DG.ViewBase;

  const csvDf1 = `seq
ATC-G-TTGC--
ATC-G-TTGC--
-TC-G-TTGC--
-TC-GCTTGC--
-TC-GCTTGC--`;


  before(async () => {
    tvList = [];
    dfList = [];
    // currentView = grok.shell.v;
  });

  after(async () => {
    // Closing opened views causes the error 'Cannot read properties of null (reading 'f')'
    // dfList.forEach((df: DG.DataFrame) => { grok.shell.closeTable(df); });
    // tvList.forEach((tv: DG.TableView) => tv.close());
    // grok.shell.v = currentView;
  });

  test('allPositions', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf1);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ.MSA');

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo')) as WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const positions: PI[] = wlViewer['positions'];

    const resAllDf1: PI[] = [
      new PI(0, '1', {'A': new PMI(2), '-': new PMI(3)}),
      new PI(1, '2', {'T': new PMI(5)}),
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(3, '4', {'-': new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {'-': new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
      new PI(7, '8', {'T': new PMI(5)}),
      new PI(8, '9', {'G': new PMI(5)}),
      new PI(9, '10', {'C': new PMI(5)}),
      new PI(10, '11', {'-': new PMI(5)}),
      new PI(11, '12', {'-': new PMI(5)}),
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const key in positions[i].freq) {
        expect(positions[i].freq[key].count, resAllDf1[i].freq[key].count);
      }
    }
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
    seqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ');

    df.filter.init((i) => {
      return i > 2;
    });
    df.filter.fireChanged();
    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'shrinkEmptyTail': true})) as WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const positions: PI[] = wlViewer['positions'];

    const resAllDf1: PI[] = [
      new PI(0, '1', {'-': new PMI(3)}),
      new PI(1, '2', {'T': new PMI(3)}),
      new PI(2, '3', {'-': new PMI(3)}),
      new PI(3, '4', {'-': new PMI(3)}),
      new PI(4, '5', {'C': new PMI(3)}),
      new PI(5, '6', {'-': new PMI(2), 'C': new PMI(1)}),
      new PI(6, '7', {'G': new PMI(3)}),
      new PI(7, '8', {'T': new PMI(3)}),
      new PI(8, '9', {'-': new PMI(3)}),
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const key in positions[i].freq) {
        expect(positions[i].freq[key].count, resAllDf1[i].freq[key].count);
      }
    }
  });

  test('positions with skipEmptyPositions option', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf1);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.DNA);
    seqCol.setTag(bioTAGS.aligned, 'SEQ.MSA');

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'skipEmptyPositions': true})) as WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const resPosList: PI[] = wlViewer['positions'];

    const tgtPosList: PI[] = [
      new PI(0, '1', {'A': new PMI(2), '-': new PMI(3)}),
      new PI(1, '2', {'T': new PMI(5)}),
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {'-': new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
      new PI(7, '8', {'T': new PMI(5)}),
      new PI(8, '9', {'G': new PMI(5)}),
      new PI(9, '10', {'C': new PMI(5)})
    ];

    expect(resPosList.length, tgtPosList.length);
    for (let posI = 0; posI < resPosList.length; posI++) {
      const resPos = resPosList[posI];
      const tgtPos = tgtPosList[posI];
      expectPositionInfo(resPos, tgtPos);
    }
  });

  test('count sequences for monomer at position', async () => {
    const df: DG.DataFrame = buildDfWithSeqCol(csvDf1, NOTATION.FASTA, ALPHABET.DNA, 'SEQ.MSA');
    const seqCol: DG.Column = df.getCol('seq');

    const tv: DG.TableView = grok.shell.addTableView(df);

    const wlViewer: WebLogoViewer = (await df.plot.fromType('WebLogo', {
      startPositionName: '3',
      endPositionName: '7',
      skipEmptyPositions: true
    })) as WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const resPosList: PI[] = wlViewer['positions'];
    const tgtPosList: PI[] = [
      new PI(2, '3', {'C': new PMI(5)}),
      new PI(4, '5', {'G': new PMI(5)}),
      new PI(5, '6', {'-': new PMI(3), 'C': new PMI(2)}),
      new PI(6, '7', {'T': new PMI(5)}),
    ];

    expect(resPosList.length, tgtPosList.length);
    for (let posI = 0; posI < resPosList.length; posI++) {
      const resPos = resPosList[posI];
      const tgtPos = tgtPosList[posI];
      expectPositionInfo(resPos, tgtPos);
    }

    const splitter: SplitterFunc = wlViewer['splitter']!;
    const atPI1: PI = resPosList[1];
    const countAt1 = countForMonomerAtPosition(df, seqCol, df.filter, splitter, 'G', atPI1);
    expect(countAt1, 5);
  });
});

function expectPositionInfo(actualPos: PI, expectedPos: PI): void {
  expect(actualPos.name, expectedPos.name);
  expectArray(Object.keys(actualPos.freq), Object.keys(expectedPos.freq));
  for (const key in actualPos.freq) {
    //
    expect(actualPos.freq[key].count, expectedPos.freq[key].count);
  }
}

function buildDfWithSeqCol(csv: string, notation: NOTATION, alphabet: ALPHABET, aligned: string): DG.DataFrame {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);

  const seqCol: DG.Column = df.getCol('seq');
  seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  seqCol.setTag(DG.TAGS.UNITS, notation);
  seqCol.setTag(bioTAGS.alphabet, alphabet);
  seqCol.setTag(bioTAGS.aligned, aligned);

  return df;
}
