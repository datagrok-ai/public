import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';
import {startAnalysis} from '../widgets/peptides';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

category('Core', () => {
  let simpleTable: DG.DataFrame;
  let simpleActivityCol: DG.Column<number>;
  let simpleAlignedSeqCol: DG.Column<string>;
  let simpleScaledCol: DG.Column<number>;

  let complexTable: DG.DataFrame;
  let complexActivityCol: DG.Column<number>;
  let complexAlignedSeqCol: DG.Column<string>;
  let complexScaledCol: DG.Column<number>;
  const alignedSequenceCol = 'AlignedSequence';

  let model: PeptidesModel | null = null;

  test('Start analysis: simple', async () => {
    const simpleActivityColName = 'IC50';
    simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
    simpleActivityCol = simpleTable.getCol(simpleActivityColName);
    simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
    simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
    simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(
      simpleActivityCol, simpleAlignedSeqCol, null, simpleTable, simpleScaledCol, C.SCALING_METHODS.MINUS_LG);
    expect(model instanceof PeptidesModel, true);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

    model!.mutationCliffsSelection = {'11': ['D']};
  });

  test('Start analysis: сomplex', async () => {
    const complexActivityColName = 'Activity';
    complexTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned_2.csv'));
    complexActivityCol = complexTable.getCol(complexActivityColName);
    complexAlignedSeqCol = complexTable.getCol('MSA');
    complexAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    complexAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.UN);
    complexAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    complexAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    complexAlignedSeqCol.setTag(C.TAGS.SEPARATOR, '/');
    complexScaledCol = scaleActivity(complexActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(
      complexActivityCol, complexAlignedSeqCol, null, complexTable, complexScaledCol, C.SCALING_METHODS.MINUS_LG);
    expect(model instanceof PeptidesModel, true);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

    if (model !== null)
      model.mutationCliffsSelection = {'13': ['-']};
  });

  test('Save and load project', async () => {
    const simpleActivityColName = 'IC50';
    simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
    simpleActivityCol = simpleTable.getCol(simpleActivityColName);
    simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
    simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
    simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
    simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);

    model = await startAnalysis(
      simpleActivityCol, simpleAlignedSeqCol, null, simpleTable, simpleScaledCol, C.SCALING_METHODS.MINUS_LG);

    let v = grok.shell.getTableView('Peptides analysis');

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => v.table!.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

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
    await awaitCheck(() => typeof grok.shell.tableView('Peptides analysis') === 'undefined', 'Table never closed', 2000);

    await sp.open();
    v = grok.shell.getTableView('Peptides analysis');

    await grok.dapi.layouts.delete(sl);
    await grok.dapi.tables.delete(sti);
    await grok.dapi.projects.delete(sp);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => v.table!.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
  });

  test('Cluster stats - Benchmark HELM 5k', async () => {
    if (!DG.Test.isInBenchmark)
      return;

    const df = (await _package.files.readBinaryDataFrames('tests/aligned_5k_2.d42'))[0];
    const activityCol = df.getCol('Activity');
    const scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    const clustersCol = df.getCol('Cluster');
    const sequenceCol = df.getCol('HELM');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    const model = await startAnalysis(
      activityCol, sequenceCol, clustersCol, df, scaledActivityCol, C.SCALING_METHODS.NONE);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

    for (let i = 0; i < 5; ++i)
      DG.time('Cluster stats', () => model?.calculateClusterStatistics());
  }, {timeout: 10000});

  test('Monomer Position stats - Benchmark HELM 5k', async () => {
    if (!DG.Test.isInBenchmark)
      return;

    const df = (await _package.files.readBinaryDataFrames('tests/aligned_5k.d42'))[0];
    const activityCol = df.getCol('Activity');
    const scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    const clustersCol = df.getCol('Cluster');
    const sequenceCol = df.getCol('HELM');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    const model = await startAnalysis(
      activityCol, sequenceCol, clustersCol, df, scaledActivityCol, C.SCALING_METHODS.NONE);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

    for (let i = 0; i < 5; ++i)
      DG.time('Monomer position stats', () => model?.calculateMonomerPositionStatistics());
  }, {timeout: 10000});

  test('Analysis start - Benchmark HELM 5k', async () => {
    if (!DG.Test.isInBenchmark)
      return;

    const df = (await _package.files.readBinaryDataFrames('tests/aligned_5k.d42'))[0];
    const activityCol = df.getCol('Activity');
    const scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    const clustersCol = df.getCol('Cluster');
    const sequenceCol = df.getCol('HELM');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);

    for (let i = 0; i < 5; ++i) {
      await DG.timeAsync('Analysis start', async () => {
        const model = await startAnalysis(
          activityCol, sequenceCol, clustersCol, df, scaledActivityCol, C.SCALING_METHODS.NONE);

        // Ensure grid finished initializing to prevent Unhandled exceptions
        let accrodionInit = false;
        grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
        await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
        await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
        await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);

        if (model)
          grok.shell.closeTable(model.df);
      });
    }
  }, {timeout: 10000});
});
