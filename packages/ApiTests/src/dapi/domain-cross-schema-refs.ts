import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Cross-schema refs (GROK-20298): a domain column pointing at a Core table the platform may
// hard-delete is a SOFT ref — no FK, the engine probes the target on single-row writes,
// filters/expands/facets hop into it under the target's row security, and a deleted target
// leaves the row readable with the dangling id. Fixture: PlatesFixture's `plate.source_query`
// (`ref: Core.queries`); skips cleanly where that column is not deployed — the
// domain-handlers.ts pattern.
category('Dapi: domain cross-schema refs', () => {
  const plates = () => grok.dapi.domains.table('plates.plate');
  const unique = (prefix: string) => `${prefix}-${Date.now()}-${Math.floor(Math.random() * 1e9)}`;
  const missingId = '00000000-0000-4000-8000-000000000000';

  async function deployed(): Promise<boolean> {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'plates'))
      return false;
    const props = await grok.dapi.domains.registry.rowProperties('plates.plate');
    return props.some((p) => p.name === 'source_query');
  }

  async function createQuery(name: string): Promise<_DG.DataQuery> {
    const dc = await grok.dapi.connections.filter('name = "Datagrok"').first();
    const q = dc.query(name, 'select 1 as x');
    q.newId();
    return await grok.dapi.queries.save(q);
  }

  async function thrown(action: () => Promise<any>): Promise<any> {
    try {
      await action();
      return null;
    } catch (e) {
      return e;
    }
  }

  test('insert: a visible query id is accepted, an unknown one is refused per row with fk', async () => {
    if (!await deployed()) {
      console.log('skipped: plates.plate.source_query is not deployed');
      return;
    }
    const query = await createQuery(unique('xref query'));
    let plateId: string | null = null;
    try {
      const [ins] = await plates().insert({barcode: unique('XREF'), source_query: query.id});
      expect(ins.created, true);
      plateId = ins.id;
      expect((await plates().get(plateId)).source_query, query.id);

      const err = await thrown(() => plates().insert({barcode: unique('XREF'), source_query: missingId}));
      expect(err instanceof DG.DomainValidationError, true,
        `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(err.status, 409);
      expect(err.code, 'validation');
      const rowError = err.rows[0].errors[0];
      expect(rowError.column, 'source_query');
      expect(rowError.code, 'fk');
      expect(rowError.message, 'Referenced row does not exist');
    } finally {
      if (plateId != null)
        await plates().delete(plateId);
      await grok.dapi.queries.delete(query);
    }
  });

  test('filter, expand and facet travel the ref into Core.queries', async () => {
    if (!await deployed()) {
      console.log('skipped: plates.plate.source_query is not deployed');
      return;
    }
    const name = unique('xref query');
    const barcode = unique('XREF');
    const query = await createQuery(name);
    let plateId: string | null = null;
    try {
      const [ins] = await plates().insert({barcode: barcode, source_query: query.id});
      plateId = ins.id;

      const filtered = await plates().query({filter: `source_query.friendlyName = "${name}"`});
      expect(filtered.length, 1, `filter through the ref: ${JSON.stringify(filtered)}`);
      expect(filtered[0].id, plateId);

      const [expanded] = await plates().query({filter: `barcode = "${barcode}"`, expand: ['source_query']});
      expect(expanded['source_query.friendlyName'], name, `expand: ${JSON.stringify(expanded)}`);

      const res = await plates().facets({facets: [
        {id: 'q', kind: 'categories', column: 'source_query', search: name},
      ]});
      const cat = (res.facets['q'] as _DG.DomainFacetCategoriesResult).categories
        .find((c) => c.value === query.id);
      expect(cat != null, true, `facet on the ref column: ${JSON.stringify(res.facets['q'])}`);
      expect(cat!.display, name);
      expect(cat!.total, 1);
    } finally {
      if (plateId != null)
        await plates().delete(plateId);
      await grok.dapi.queries.delete(query);
    }
  });

  test('deleted query: the row keeps the dangling id, expand yields nothing, filter finds no rows', async () => {
    if (!await deployed()) {
      console.log('skipped: plates.plate.source_query is not deployed');
      return;
    }
    const name = unique('xref query');
    const barcode = unique('XREF');
    const query = await createQuery(name);
    let plateId: string | null = null;
    try {
      const [ins] = await plates().insert({barcode: barcode, source_query: query.id});
      plateId = ins.id;
      await grok.dapi.queries.delete(query);

      expect((await plates().get(plateId)).source_query, query.id, 'the dangling id must stay on the row');
      const [expanded] = await plates().query({filter: `barcode = "${barcode}"`, expand: ['source_query']});
      expect(expanded['source_query.friendlyName'] == null, true, `expand: ${JSON.stringify(expanded)}`);
      expect((await plates().query({filter: `source_query.friendlyName = "${name}"`})).length, 0);
    } finally {
      if (plateId != null)
        await plates().delete(plateId);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
