import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {after, before, category, test, expect, expectObject, delay} from '@datagrok-libraries/utils/src/test';

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
    currentView = grok.shell.v;
  });

  after(async () => {
    dfList.forEach((df: DG.DataFrame) => { grok.shell.closeTable(df);});
    tvList.forEach((tv: DG.TableView) => tv.close());
    grok.shell.v = currentView;
  });

  test('allPositions', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf1);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const seqCol: DG.Column = df.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.setTag(DG.TAGS.UNITS, bio.NOTATION.FASTA);
    seqCol.setTag(bio.TAGS.alphabet, bio.ALPHABET.DNA);
    seqCol.setTag(bio.TAGS.aligned, 'SEQ.MSA');

    const wlViewer: bio.WebLogoViewer = (await df.plot.fromType('WebLogo')) as bio.WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const positions: bio.PositionInfo[] = wlViewer['positions'];

    const resAllDf1: bio.PositionInfo[] = [
      new bio.PositionInfo('1', {'A': new bio.PositionMonomerInfo(2), '-': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('2', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('3', {'C': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('4', {'-': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('5', {'G': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('6', {'-': new bio.PositionMonomerInfo(3), 'C': new bio.PositionMonomerInfo(2)}),
      new bio.PositionInfo('7', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('8', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('9', {'G': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('10', {'C': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('11', {'-': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('12', {'-': new bio.PositionMonomerInfo(5)}),
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const key in positions[i].freq) {
        expect(positions[i].freq[key].count, resAllDf1[i].freq[key].count);
      }
    }
  });

  test('positions with shrinkEmptyTail option true (filterd)', async () => {
    let csvDf2 = `seq 
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
    seqCol.setTag(DG.TAGS.UNITS, bio.NOTATION.FASTA);
    seqCol.setTag(bio.TAGS.alphabet, bio.ALPHABET.DNA);
    seqCol.setTag(bio.TAGS.aligned, 'SEQ');

    df.filter.init((i) => {
      return i > 2;
    });
    df.filter.fireChanged();
    const wlViewer: bio.WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'shrinkEmptyTail': true})) as bio.WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const positions: bio.PositionInfo[] = wlViewer['positions'];

    const resAllDf1: bio.PositionInfo[] = [
      new bio.PositionInfo('1', {'-': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('2', {'T': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('3', {'-': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('4', {'-': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('5', {'C': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('6', {'-': new bio.PositionMonomerInfo(2), 'C': new bio.PositionMonomerInfo(1)}),
      new bio.PositionInfo('7', {'G': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('8', {'T': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('9', {'-': new bio.PositionMonomerInfo(3)}),
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
    seqCol.setTag(DG.TAGS.UNITS, bio.NOTATION.FASTA);
    seqCol.setTag(bio.TAGS.alphabet, bio.ALPHABET.DNA);
    seqCol.setTag(bio.TAGS.aligned, 'SEQ.MSA');

    const wlViewer: bio.WebLogoViewer = (await df.plot.fromType('WebLogo',
      {'skipEmptyPositions': true})) as bio.WebLogoViewer;
    tv.dockManager.dock(wlViewer.root, DG.DOCK_TYPE.DOWN);

    tvList.push(tv);
    dfList.push(df);

    const positions: bio.PositionInfo[] = wlViewer['positions'];

    const resAllDf1: bio.PositionInfo[] = [
      new bio.PositionInfo('1', {'A': new bio.PositionMonomerInfo(2), '-': new bio.PositionMonomerInfo(3)}),
      new bio.PositionInfo('2', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('3', {'C': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('5', {'G': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('6', {'-': new bio.PositionMonomerInfo(3), 'C': new bio.PositionMonomerInfo(2)}),
      new bio.PositionInfo('7', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('8', {'T': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('9', {'G': new bio.PositionMonomerInfo(5)}),
      new bio.PositionInfo('10', {'C': new bio.PositionMonomerInfo(5)})
    ];

    expect(positions.length, resAllDf1.length);

    for (let i = 0; i < positions.length; i++) {
      expect(positions[i].name, resAllDf1[i].name);
      for (const key in positions[i].freq) {
        expect(positions[i].freq[key].count, resAllDf1[i].freq[key].count);
      }
    }
  });
});
