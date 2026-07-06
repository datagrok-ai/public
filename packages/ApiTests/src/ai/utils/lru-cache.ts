import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

// Tests DG.LruCache.
category('AI: Utils: LRU cache', () => {
  test('set/get/has round-trip and miss returns undefined', async () => {
    const c = new DG.LruCache<string, number>();
    expect(c.has('a'), false);
    expect(c.get('a') === undefined, true);
    c.set('a', 1);
    c.set('b', 2);
    expect(c.has('a'), true);
    expect(c.has('b'), true);
    expect(c.get('a'), 1);
    expect(c.get('b'), 2);
    expect(c.has('missing'), false);
    expect(c.get('missing') === undefined, true);
  });

  test('set on existing key updates value without growing capacity', async () => {
    const c = new DG.LruCache<string, number>(2);
    c.set('a', 1);
    c.set('a', 2);
    c.set('a', 3);
    expect(c.get('a'), 3);
    c.set('b', 10);
    expect(c.has('a'), true);
    expect(c.has('b'), true);
    expect(c.get('b'), 10);
  });

  test('eviction: oldest key dropped when capacity exceeded', async () => {
    const c = new DG.LruCache<string, number>(2);
    c.set('a', 1);
    c.set('b', 2);
    c.set('c', 3);
    expect(c.has('a'), false);
    expect(c.has('b'), true);
    expect(c.has('c'), true);
    expect(c.get('a') === undefined, true);
    expect(c.get('b'), 2);
    expect(c.get('c'), 3);
  });

  test('onItemEvicted callback fires with the evicted value', async () => {
    const c = new DG.LruCache<string, string>(2);
    const evicted: string[] = [];
    c.onItemEvicted = (v: string) => evicted.push(v);
    c.set('a', 'va');
    c.set('b', 'vb');
    expectArray(evicted, []);
    c.set('c', 'vc');
    expectArray(evicted, ['va']);
    c.set('d', 'vd');
    expectArray(evicted, ['va', 'vb']);
  });

  test('LRU ordering: get bumps key so it is not evicted next', async () => {
    const c = new DG.LruCache<string, number>(2);
    c.set('a', 1);
    c.set('b', 2);
    expect(c.get('a'), 1);
    c.set('c', 3);
    expect(c.has('a'), true);
    expect(c.has('b'), false);
    expect(c.has('c'), true);
  });

  test('LRU ordering: re-set bumps key so it is not evicted next', async () => {
    const c = new DG.LruCache<string, number>(2);
    c.set('a', 1);
    c.set('b', 2);
    c.set('a', 11);
    c.set('c', 3);
    expect(c.has('a'), true);
    expect(c.has('b'), false);
    expect(c.has('c'), true);
    expect(c.get('a'), 11);
  });

  test('getOrCreate returns cached on hit and invokes factory on miss', async () => {
    const c = new DG.LruCache<string, number>(4);
    let calls = 0;
    const v1 = c.getOrCreate('k', (k) => {
      calls++;
      return k.length + 100;
    });
    expect(v1, 101);
    expect(calls, 1);
    const v2 = c.getOrCreate('k', (_k) => {
      calls++;
      return 999;
    });
    expect(v2, 101);
    expect(calls, 1);
    expect(c.has('k'), true);
    expect(c.get('k'), 101);
  });
}, {owner: 'agolovko@datagrok.ai'});
