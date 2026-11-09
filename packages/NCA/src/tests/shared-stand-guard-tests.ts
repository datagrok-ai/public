import {category, test, expect, tests} from '@datagrok-libraries/test/src/test';

import {
  STORE_WRITING_TESTS, SHARED_STAND_SKIP_REASON, SHARED_STAND_WRITE_REFUSED, isLoopbackServer, isStoreWriting,
  skipStoreWritingTestsOnSharedStand, refuseStoreWrites,
} from './shared-stand-guard';
import {MemoryNcaRepository} from './memory-nca-repository';

/** The registry as names only — the guard is exercised on a copy, never on the live registry. */
function registryCopy(): Record<string, {tests: {name: string; options?: {skipReason?: string}}[]}> {
  const out: Record<string, {tests: {name: string; options?: {skipReason?: string}}[]}> = {};
  for (const [cat, entry] of Object.entries(tests as Record<string, {tests?: {name: string}[]}>))
    out[cat] = {tests: (entry.tests ?? []).map((t) => ({name: t.name}))};
  return out;
}

category('Tests: shared-stand guard', () => {
  test('loopback servers run the store tests; any other host does not', async () => {
    expect(isLoopbackServer('http://localhost:8080/api'), true, 'local stand');
    expect(isLoopbackServer('http://127.0.0.1:8080/api'), true, 'CI stand (packages.yaml apiUrl)');
    expect(isLoopbackServer('https://dev.datagrok.ai/api'), false, 'dev');
    expect(isLoopbackServer('https://localhost.example.com/api'), false, 'a name that only starts with localhost');
    expect(isLoopbackServer('not a url'), false, 'unparseable → not loopback (skip, the safe side)');
  });

  test('every rule matches a registered test (a renamed category cannot slip past)', async () => {
    const reg = registryCopy();
    for (const rule of STORE_WRITING_TESTS) {
      const hit = Object.entries(reg).some(([cat, e]) =>
        rule.category.test(cat) && e.tests.some((t) => rule.test === undefined || rule.test.test(t.name)));
      expect(hit, true, `rule ${rule.category}${rule.test ? ' / ' + rule.test : ''} matches no registered test`);
    }
  });

  test('on a shared server the store tests are skipped with the reason; the rest are untouched', async () => {
    const reg = registryCopy();
    const marked = skipStoreWritingTestsOnSharedStand(reg, 'https://dev.datagrok.ai/api');
    let expected = 0;
    for (const [cat, e] of Object.entries(reg)) {
      for (const t of e.tests) {
        const writes = isStoreWriting(cat, t.name);
        if (writes) expected++;
        expect(t.options?.skipReason === SHARED_STAND_SKIP_REASON, writes, `${cat} / ${t.name}`);
      }
    }
    expect(marked, expected, 'marked count');
    expect(marked > 0, true, 'some store tests exist');
    expect(isStoreWriting('Tests: shared-stand guard', 'x'), false, 'the guard tests themselves run everywhere');
    expect(skipStoreWritingTestsOnSharedStand(registryCopy(), 'http://localhost:8080/api'), 0, 'loopback marks none');
    return `${marked} store-writing tests`;
  });

  test('the backstop repository refuses writes and still reads', async () => {
    const inner = new MemoryNcaRepository();
    const guarded = refuseStoreWrites(inner);
    let error = '';
    try {
      await guarded.createStudy('zz-guard');
    } catch (e) {
      error = e instanceof Error ? e.message : String(e);
    }
    expect(error.startsWith(SHARED_STAND_WRITE_REFUSED), true, `write refused, got "${error}"`);
    expect((await inner.listStudies()).length, 0, 'nothing reached the store');
    expect(Array.isArray(await guarded.listStudies()), true, 'a read passes through');
  });
});
