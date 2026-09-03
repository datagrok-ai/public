import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {pltsDb} from '../generated/db';
import {initPlates, createPlateTemplate, deletePlateTemplate, getOrCreateProperty,
  allProperties, plateTemplates, plateTypes, savePlate, getPlateById,
  createAnalysisRun, saveAnalysisResult} from '../plates/plates-crud';
import {Plate} from '../plate/plate';

// Persistence tests for the plts domain schema (databases/plts/schema.json).
// Requires a stand where Plates:SetupPltsSchema has run (grants + dictionary seeds).
category('Plates: dictionaries CRUD', () => {
  const unique = (prefix: string) => `${prefix}-${Date.now() % 1e10}-${Math.floor(Math.random() * 1e3)}`;

  test('initPlates loads seeded dictionaries', async () => {
    await initPlates(true);
    const type96 = plateTypes.find((pt) => pt.rows === 8 && pt.cols === 12);
    expect(type96 != null, true, 'seeded 96-well plate type not found');
    expect(type96!.id !== '', true, 'plate type id not loaded');
    const wellRole = allProperties.find((p) => p.name === 'Well Role' && p.scope === 'well');
    expect(wellRole != null, true, 'seeded Well Role property not found');
    expect(JSON.parse(wellRole!.choices as any as string).includes('DMSO'), true);
  });

  test('property get-or-create dedups by name and scope', async () => {
    const name = unique('prop');
    const first = await getOrCreateProperty(name, DG.TYPE.FLOAT, 'well');
    try {
      const second = await getOrCreateProperty(name, DG.TYPE.FLOAT, 'well');
      expect(second.id, first.id, 'same name+scope must resolve to the same property');
      const [report] = await pltsDb.propertieses.insert({name: name, type: 'double', scope: 'well'});
      expect(report.created, false, 'business key must dedup a direct duplicate insert');
      expect(report.existingId, first.id);
    } finally {
      await pltsDb.propertieses.delete(first.id);
      await initPlates(true);
    }
  });

  test('template lifecycle: create, reload, delete', async () => {
    const templateName = unique('template');
    const plateProp = unique('plate-prop');
    const template = await createPlateTemplate({
      name: templateName,
      description: 'test template',
      plateProperties: [{name: plateProp, type: DG.TYPE.STRING, is_required: true, default_value: 'x'}],
      wellProperties: [{name: 'Concentration', type: DG.TYPE.FLOAT}],
    });
    try {
      expect(template.id !== '', true, 'template id not assigned');
      expect(template.required_props.length, 1);

      await initPlates(true);
      const reloaded = plateTemplates.find((t) => t.id === template.id)!;
      expect(reloaded != null, true, 'template not found after reload');
      expect(reloaded.plateProperties.length, 1);
      expect(reloaded.plateProperties[0].name, plateProp);
      expect((reloaded.plateProperties[0] as any).is_required, true);
      expect((reloaded.plateProperties[0] as any).default_value, 'x');
      const concentration = allProperties.find((p) => p.name === 'Concentration' && p.scope === 'well')!;
      expect(reloaded.wellProperties[0].id, concentration.id, 'seeded property must be reused, not duplicated');
    } finally {
      await deletePlateTemplate(template);
      expect(await pltsDb.templatePropertieses.count(DG.cond('template_id', '=', template.id)), 0,
        'template_properties rows must cascade with the template');
      const prop = allProperties.find((p) => p.name === plateProp && p.scope === 'plate');
      if (prop)
        await pltsDb.propertieses.delete(prop.id);
      await initPlates(true);
    }
  });
});

category('Plates: persistence', () => {
  const unique = (prefix: string) => `${prefix}-${Date.now() % 1e10}-${Math.floor(Math.random() * 1e3)}`;

  test('plate save/read round-trip', async () => {
    await initPlates(true);
    const template = await createPlateTemplate({
      name: unique('rt-template'),
      plateProperties: [{name: unique('rt-chemist'), type: DG.TYPE.STRING}],
      wellProperties: [{name: 'Concentration', type: DG.TYPE.FLOAT}],
    });
    const chemistProp = template.plateProperties[0].name!;
    const plate = new Plate(8, 12);
    plate.details = {[chemistProp]: 'John Doe'};
    plate.data.columns.addNew('Concentration', DG.COLUMN_TYPE.FLOAT).init((i) => i);
    plate.plateTemplateId = template.id;

    await savePlate(plate);
    try {
      expect(plate.id != null && plate.id !== '', true, 'plate id not assigned');
      const saved = await pltsDb.plateses.get(plate.id!);
      expect(saved.template_id, template.id, 'template link not persisted');
      expect(await pltsDb.plateDetailses.count(DG.cond('plate_id', '=', plate.id!)), 1);

      const reloaded = await getPlateById(plate.id!);
      expect(reloaded.rows, 8);
      expect(reloaded.cols, 12);
      expect(reloaded.data.col('Concentration')!.get(50), 50);

      const runId = await createAnalysisRun(plate.id!, 'rt-analysis', ['GRK-1']);
      const ic50 = await getOrCreateProperty(unique('rt-ic50'), DG.TYPE.FLOAT);
      await saveAnalysisResult({runId: runId, propertyId: ic50.id, propertyName: ic50.name,
        propertyType: DG.TYPE.FLOAT, value: 42.5, groupCombination: ['GRK-1']});
      expect(await pltsDb.analysisResultses.count(DG.cond('analysis_run_id', '=', runId)), 1);
      await pltsDb.analysisRunses.delete(runId);
      expect(await pltsDb.analysisResultses.count(DG.cond('analysis_run_id', '=', runId)), 0,
        'analysis results must cascade with the run');
      await pltsDb.propertieses.delete(ic50.id);
    } finally {
      if (plate.id)
        await pltsDb.plateses.delete(plate.id);
      await deletePlateTemplate(template);
      const prop = allProperties.find((p) => p.name === chemistProp && p.scope === 'plate');
      if (prop)
        await pltsDb.propertieses.delete(prop.id);
      await initPlates(true);
    }
  });
});

// Write-throughput measurements: the sequential test reproduces the legacy
// one-INSERT-per-value pattern, the batch test is the migrated path — the two
// durations side by side are the speed-up. Run with `grok test --benchmark`
// for the full sizes; normal runs use small ones.
category('Plates: benchmarks', () => {
  // Walks a 32×48 grid — up to 1536 rows with distinct well positions (a full 1536-well plate).
  const wellValueRows = (plateId: string, propertyId: string, count: number) =>
    Array.from({length: count}, (_, i) => ({
      plate_id: plateId, row: Math.floor(i / 48), col: i % 48,
      property_id: propertyId, value_num: Math.random() * 100,
    }));

  const withPlate = async (action: (plateId: string, propertyId: string) => Promise<void>) => {
    await initPlates(true);
    const type96 = plateTypes.find((pt) => pt.rows === 8 && pt.cols === 12)!;
    const property = await getOrCreateProperty(`bench-${Date.now()}`, DG.TYPE.FLOAT, 'well');
    const [plate] = await pltsDb.plateses.insert({plate_type_id: type96.id, barcode: `bench-${Date.now()}`});
    try {
      await action(plate.id, property.id);
    } finally {
      await pltsDb.plateses.delete(plate.id);
      await pltsDb.propertieses.delete(property.id);
      await initPlates(true);
    }
  };

  test('well values: single-row inserts (legacy write pattern)', async () => {
    const count = DG.Test.isInBenchmark ? 192 : 24;
    await withPlate(async (plateId, propertyId) => {
      for (const row of wellValueRows(plateId, propertyId, count))
        await pltsDb.plateWellValueses.insert(row);
    });
  }, {benchmark: true, benchmarkTimeout: 120000});

  test('well values: one batch (migrated write pattern)', async () => {
    const count = DG.Test.isInBenchmark ? 1536 : 96;
    await withPlate(async (plateId, propertyId) => {
      const report = await pltsDb.plateWellValueses.batch(wellValueRows(plateId, propertyId, count));
      expect(report.inserted, count);
    });
  }, {benchmark: true, benchmarkTimeout: 120000});

  test('initPlates: cold dictionary load', async () => {
    await initPlates(true);
  }, {benchmark: true});
});
