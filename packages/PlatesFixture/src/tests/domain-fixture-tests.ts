import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Engine-correctness fixture for the 'plates' domain schema declared in
// databases/plates/schema.json: all three security modes, constraints, audit.
category('Domains: plates fixture', () => {
  const types = () => grok.dapi.domains.table('plates.plate_type');
  const plates = () => grok.dapi.domains.table('plates.plate');
  const wells = () => grok.dapi.domains.table('plates.plate_well');
  const readings = () => grok.dapi.domains.table('plates.well_reading');
  const unique = (prefix: string) => `${prefix}-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  const throws = async (action: () => Promise<any>): Promise<string> => {
    try {
      await action();
    } catch (e: any) {
      return `${e.message ?? e}`;
    }
    return '';
  };

  test('schema registration', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    const s = schemas.find((x) => x.name === 'plates');
    expect(s != null, true, 'plates schema not registered');
    const modes = new Map(s!.tables.map((t) => [t.name, t.securityMode]));
    expect(modes.get('plate_type'), 'table');
    expect(modes.get('plate'), 'row');
    expect(modes.get('plate_well'), 'master');
    expect(modes.get('well_reading'), 'master');
  });

  test('table mode: crud, audit off', async () => {
    const [ins] = await types().insert({name: unique('type'), rows: 8, cols: 12});
    expect(ins.created, true);
    try {
      await types().update(ins.id, {rows: 16}, {version: 1});
      const row = await types().get(ins.id);
      expect(row.rows, 16);
      const audit = await types().audit(ins.id);
      expect(audit.length, 0, 'audit-off table must write no audit rows');
    } finally {
      await types().delete(ins.id);
    }
  });

  test('row mode: jsonb schemas, dedup, audit, promote', async () => {
    const barcode = unique('BC');
    const [ins] = await plates().insert(
      {barcode: barcode, row_count: 8, col_count: 12, organism: 'human', price: 9.99});
    try {
      const [row] = await plates().query({filter: `barcode = "${barcode}"`});
      expect(row.organism, 'human', 'jsonb property schema key not returned');
      expect(row.price, 9.99);

      const [dup] = await plates().insert({barcode: barcode});
      expect(dup.status, 'duplicate');
      expect(dup.existingId, ins.id);

      await plates().update(ins.id, {organism: 'mouse'}, {version: 1});
      const audit = await plates().audit(ins.id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].after.data.organism, 'mouse', 'jsonb keys must appear in the audit diff');

      const promoted = await plates().promote(ins.id);
      expect(promoted.promoted, true);
    } finally {
      await plates().delete(ins.id);
    }
    expect(await plates().get(ins.id), null, 'soft-deleted plate still visible');
  });

  test('master mode: wells cascade with their plate', async () => {
    const [plate] = await plates().insert({barcode: unique('BC')});
    await wells().insert([
      {plate_id: plate.id, row: 0, col: 0, volume: 10.5, sample_state: 'filled'},
      {plate_id: plate.id, row: 0, col: 1},
    ]);
    expect((await wells().query({filter: `plate_id = "${plate.id}"`})).length, 2);
    await plates().delete(plate.id);
    expect((await wells().query({filter: `plate_id = "${plate.id}"`})).length, 0,
      'cascade soft-delete did not remove wells');
  });

  test('depth-2 master mode: readings follow the plate through the chain', async () => {
    const [plate] = await plates().insert({barcode: unique('BC')});
    const [well] = await wells().insert({plate_id: plate.id, row: 0, col: 0});
    await readings().insert([
      {well_id: well.id, wavelength: 450, value: 0.42},
      {well_id: well.id, wavelength: 620, value: 0.13},
    ]);
    expect((await readings().query({filter: `well_id = "${well.id}"`})).length, 2);
    await plates().delete(plate.id);
    expect((await readings().query({filter: `well_id = "${well.id}"`})).length, 0,
      'cascade soft-delete did not propagate to depth-2 readings');
  });

  test('constraints: required, min, choices', async () => {
    expect((await throws(() => plates().insert({row_count: 8}))) !== '', true,
      'missing required barcode accepted');
    expect((await throws(() => plates().insert({barcode: unique('BC'), row_count: 0}))) !== '', true,
      'row_count below min accepted');
    const [plate] = await plates().insert({barcode: unique('BC')});
    try {
      expect((await throws(() => wells().insert({plate_id: plate.id, row: 0, col: 0, volume: -1}))) !== '',
        true, 'negative volume accepted');
      expect((await throws(() => wells().insert(
        {plate_id: plate.id, row: 0, col: 0, sample_state: 'bogus'}))) !== '', true,
      'value outside choices accepted');
    } finally {
      await plates().delete(plate.id);
    }
  });

  test('referential restrict on plate_type', async () => {
    const [type] = await types().insert({name: unique('type')});
    const [plate] = await plates().insert({barcode: unique('BC'), plate_type_id: type.id});
    const error = await throws(() => types().delete(type.id));
    expect(error !== '', true, 'restrict delete of a referenced plate_type succeeded');
    await plates().delete(plate.id);
    await types().delete(type.id);
    expect(await types().get(type.id), null);
  });
}, {owner: 'askalkin@datagrok.ai'});
