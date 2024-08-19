import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';
import {startAnalysis} from '../widgets/peptides';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {MonomerPosition} from '../viewers/sar-viewer';

category('Core', () => {
  const alignedSequenceCol = 'AlignedSequence';

  let model: PeptidesModel | null = null;
  test('Start analysis: simple', async () => {
    const simpleActivityColName = 'IC50';
    const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
    const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
    const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
    simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
    simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    const simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(
      simpleActivityCol, simpleAlignedSeqCol, null, simpleTable, simpleScaledCol, C.SCALING_METHODS.MINUS_LG);
    expect(model instanceof PeptidesModel, true, 'Model is null');

    const MPViewer = model!.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition | null;
    MPViewer!.mutationCliffsSelection = {'11': ['D']};
    await delay(3000);
  });

  test('Start analysis: сomplex', async () => {
    const complexActivityColName = 'Activity';
    const complexTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned_2.csv'));
    const complexActivityCol = complexTable.getCol(complexActivityColName);
    const complexAlignedSeqCol = complexTable.getCol('MSA');
    complexAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    complexAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.UN);
    complexAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    complexAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    complexAlignedSeqCol.setTag(C.TAGS.SEPARATOR, '/');
    const complexScaledCol = scaleActivity(complexActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(
      complexActivityCol, complexAlignedSeqCol, null, complexTable, complexScaledCol, C.SCALING_METHODS.MINUS_LG);
    expect(model instanceof PeptidesModel, true, 'Model is null');

    const MPViewer = model!.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition | null;
    MPViewer!.mutationCliffsSelection = {'13': ['-']};
    await delay(3000);
  });

  test('Save and load project', async () => {
    const simpleActivityColName = 'IC50';
    const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
    const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
    const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
    simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
    simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    const simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(simpleActivityCol, simpleAlignedSeqCol, null, simpleTable, simpleScaledCol,
      C.SCALING_METHODS.MINUS_LG);

    let v = grok.shell.getTableView('Peptides analysis');

    const d = v.dataFrame;
    const layout = v.saveLayout();
    const tableInfo = d.getTableInfo();

    const project = DG.Project.create();
    project.name = 'Peptides project unique test';
    project.addChild(tableInfo);
    project.addChild(layout);
    const sl = await grok.dapi.layouts.save(layout);
    await grok.dapi.tables.uploadDataFrame(d);
    const sti = await grok.dapi.tables.save(tableInfo);
    const sp = await grok.dapi.projects.save(project);

    v.close();
    await awaitCheck(() => typeof grok.shell.tableView('Peptides analysis') === 'undefined',
      'Table never closed', 3000);

    await sp.open();
    v = grok.shell.getTableView('Peptides analysis');

    await grok.dapi.layouts.delete(sl);
    await grok.dapi.tables.delete(sti);
    await grok.dapi.projects.delete(sp);
  }, {skipReason: 'ViewLayout should become ViewInfo in 1.18.'});
}, {clear: true});
