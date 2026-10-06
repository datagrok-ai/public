/* The platform seam the sources run on (WO-14), against the dg stub rather than fakes: what the
   fake `funcRunner` in sources.test.js cannot see is everything this file does with a real-shaped
   `DG.Func`. The design-time row cap lives here — `TableQuery.limit` is read while the call
   EXECUTES, so a cap restored after `prepare` is no cap at all. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

const DG = await import('datagrok-api/dg');
await import('../src/dg/shell/source-backends.js');
const {backends} = await import('../src/sources/backends.js');

/** A query that reports the cap it was actually executed under, after an await. */
function query(Kind, rows = 830) {
  const seen = [];
  const func = new Kind('Orders', {
    outputs: [new Property('result', 'dataframe')],
    run: async () => {
      await new Promise((resolve) => setTimeout(resolve, 5));
      seen.push(func.limit);
      return {result: {rowCount: func.limit ?? rows}};
    },
  });
  DG.Func.registry = [func];
  return {func, seen};
}

test('the design-time cap is still set while the call runs, and restored after it', async () => {
  const {func, seen} = query(DG.TableQuery);
  const out = await backends.funcRunner.run('Orders', {}, undefined, {maxRows: 100});
  assert.deepEqual(seen, [100], 'the limit the call executed under');
  assert.equal(out.result.rowCount, 100);
  assert.equal(func.limit, undefined, 'the shared entity keeps no design-time state');
});

test('a run with no cap leaves the query alone', async () => {
  const {seen} = query(DG.TableQuery);
  const out = await backends.funcRunner.run('Orders', {}, undefined, undefined);
  assert.deepEqual(seen, [undefined]);
  assert.equal(out.result.rowCount, 830);
});

test('a plain query has no limit to set — it runs uncapped (OR-1)', async () => {
  const {seen} = query(DG.DataQuery);
  const out = await backends.funcRunner.run('Orders', {}, undefined, {maxRows: 100});
  assert.deepEqual(seen, [undefined]);
  assert.equal(out.result.rowCount, 830, 'the recorded gap: no FuncCall-level cap exists');
});

test('two overlapping capped runs leave the shared query at its own limit', async () => {
  const {func, seen} = query(DG.TableQuery);
  func.limit = 500;
  const first = backends.funcRunner.run('Orders', {}, undefined, {maxRows: 100});
  const second = backends.funcRunner.run('Orders', {}, undefined, {maxRows: 100});
  await Promise.all([first, second]);
  assert.deepEqual(seen, [100, 100], 'neither run was uncapped by the other finishing');
  assert.equal(func.limit, 500,
    'the limit lives on the shared entity: a design-time run must not outlive itself');
});

test('find reports the kind the design-time policy defaults on', () => {
  query(DG.TableQuery);
  assert.equal(backends.funcRunner.find('Orders').kind, 'table-query',
    'a visual query — the only kind the cap and the design-time auto-run reach');
  query(DG.DataQuery);
  assert.equal(backends.funcRunner.find('Orders').kind, 'query');
  query(DG.Script);
  assert.equal(backends.funcRunner.find('Orders').kind, 'script');
  query(DG.Func);
  assert.equal(backends.funcRunner.find('Orders').kind, 'function');
  assert.equal(backends.funcRunner.find('Nope'), null);
});
