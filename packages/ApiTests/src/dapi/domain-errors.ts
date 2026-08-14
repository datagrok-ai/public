import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// WO-1 (GROK-20598): the typed domain error surface — server failures reach JS
// as DomainError subclasses discriminated on the structured body's error code,
// never by message matching. Fixtures: the package's own 'apitests' schema
// (databases/apitests/schema.json); the restrict probe rides the Grit schema
// (issue.project_id is the only deployed restrict ref) and skips cleanly where
// Grit is not deployed — the domain-handlers.ts pattern.
category('Dapi: domain errors', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const sku = () => `ERR-${Date.now()}-${Math.floor(Math.random() * 1e9)}`;

  async function thrown(action: () => Promise<any>): Promise<any> {
    try {
      await action();
      return null;
    } catch (e) {
      return e;
    }
  }

  test('version conflict: DomainVersionConflictError carries both versions', async () => {
    const [ins] = await items().insert({sku: sku(), name: 'vc probe', quantity: 1});
    try {
      await items().update(ins.id, {quantity: 2}, {version: 1}); // → version 2
      const err = await thrown(() => items().update(ins.id, {quantity: 3}, {version: 1}));
      expect(err instanceof DG.DomainVersionConflictError, true,
        `expected DomainVersionConflictError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(err.code, 'version-conflict');
      expect(err.currentVersion, 2);
      expect(err.expectedVersion, 1);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('restrict delete: DomainRestrictError', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'grit')) {
      console.log('skipped: grit schema not deployed');
      return;
    }
    const projects = grok.dapi.domains.table('grit.project');
    const issues = grok.dapi.domains.table('grit.issue');
    const [project] = await projects.insert(
      {key: `ERR${Date.now() % 10000000}`, name: 'Restrict probe'});
    const [issue] = await issues.insert(
      {project_id: project.id, number: 1, title: 'blocks project deletion'});
    try {
      const err = await thrown(() => projects.delete(project.id));
      expect(err instanceof DG.DomainRestrictError, true,
        `expected DomainRestrictError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(err.code, 'restrict');
    } finally {
      await issues.delete(issue.id);
      await projects.delete(project.id);
    }
  });

  test('malformed filter: DomainFilterError with status 400', async () => {
    const err = await thrown(() => items().query({filter: 'no_such_column = "x"'}));
    expect(err instanceof DG.DomainFilterError, true,
      `expected DomainFilterError, got ${err?.constructor?.name}: ${err?.message}`);
    expect(err.status, 400);
  });

  test('errorOnDuplicate insert: DomainValidationError.isDuplicate', async () => {
    const s = sku();
    const [first] = await items().insert({sku: s, name: 'dup probe'});
    try {
      const err = await thrown(() =>
        items().insert({sku: s, name: 'dup probe 2'}, {errorOnDuplicate: true}));
      expect(err instanceof DG.DomainValidationError, true,
        `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(err.isDuplicate, true, JSON.stringify(err?.body));
    } finally {
      await items().delete(first.id);
    }
  });
});
