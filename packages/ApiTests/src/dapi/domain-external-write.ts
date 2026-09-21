import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// The typed client over a WRITABLE external binding: the `extwlive` fixture
// (core/docs/features/ems/external-bindings/fixtures/extwlive), which the tester publishes
// before a run. Rows are addressed by their business key (`tid`), and an edit is guarded by
// `expected` — the columns' values as last read — never by `version`. The category skips
// cleanly while the fixture is not registered (a development convenience only: the release
// evidence runs it positively).
category('Dapi: domain external write', () => {
  const things = () => grok.dapi.domains.table('extwlive.thing');
  let skip: string | null = null;

  before(async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'extwlive'))
      skip = 'the extwlive fixture is not registered';
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  const isA = (err: any, ctor: any) =>
    expect(err instanceof ctor, true, `expected ${ctor.name}, got ${err?.constructor?.name}: ${err?.message}`);

  test('update with expected: a matching guard writes, a stale one names what moved', async () => {
    if (skipped())
      return;
    const tid = 900000 + (Date.now() % 90000);
    const [ins] = await things().insert({tid, note: `guard ${tid}`, n: 1});
    try {
      expect(ins.id, `${tid}`, 'an external row is addressed by the key the caller supplied');
      await things().update(ins.id, {n: 2}, {expected: {n: 1}});
      expect((await things().get(ins.id)).n, 2);

      let err = await thrown(() => things().update(ins.id, {n: 3}, {expected: {n: 1}}));
      isA(err, DG.DomainVersionConflictError);
      expect(err.code, 'version-conflict');
      expect(err.id, ins.id);
      expect(err.body.expected.n, 1, JSON.stringify(err.body));
      expect(err.body.current.n, 2, JSON.stringify(err.body));

      err = await thrown(() => things().update(ins.id, {n: 4}, {expected: {price: 1.5}}));
      isA(err, DG.DomainUnsupportedError);
      expect(err.op, 'expected', 'a float cannot guard');

      err = await thrown(() => things().update(ins.id, {n: 4}, {version: 1}));
      isA(err, DG.DomainValidationError);
      expect(err.message.includes('applies to platform-stored tables'), true, err.message);
      expect((await things().get(ins.id)).n, 2, 'a refused guard writes nothing');

      await things().delete(ins.id);
      err = await thrown(() => things().update(ins.id, {n: 5}, {expected: {n: 2}}));
      isA(err, DG.DomainNotFoundError);
    } finally {
      try {
        await things().delete(ins.id);
      } catch (e) {
        if (!(e instanceof DG.DomainNotFoundError))
          throw e;
      }
    }
  });
});
