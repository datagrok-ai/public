/**
 * Keep the store-writing tests off shared servers.
 *
 * These tests write real `nca` rows and register probe users. A delete in EMS is soft (`is_deleted`)
 * and a user cannot be deleted at all, so every run leaves rows and blocked users behind (measured on
 * the local stand 2026-11-03: 1189 soft-deleted studies, 38 probe users). On the local stand and in
 * CI (a throw-away stand per job, `http://127.0.0.1:8080`) that is harmless; on dev or any other
 * shared server it is not. So: run them only against a loopback server, skip them elsewhere with a
 * reason — the same tests, the same code, no in-memory substitute that would pass for the wrong reason.
 */

/** A test the guard may skip: a category, and optionally the test names within it. */
interface StoreWritingRule {
  readonly category: RegExp;
  /** Only these tests of the category; the whole category when absent. */
  readonly test?: RegExp;
}

/** The tests that write to the live store (EMS rows, DomainFiles snapshots, probe users). Matched by
 *  reading each test's store use (2026-11-03); a new store-writing test must be added here. */
export const STORE_WRITING_TESTS: readonly StoreWritingRule[] = [
  {category: /^NCA EMS: (platform|persisted goldens)$/},
  {category: /^NCA storage contract: EMS$/},
  {category: /^NCA storage: (EMS layout|request counts|permissions)/},
  {category: /^Pipeline functions: (save to NCA Studio|source-replaced guard)/},
  {category: /^Workflow stages: (router \+ hub|standalone store)/},
  {category: /^Workflow stages: live view/, test: /^BUG-16:/},
  {category: /^Adapter: sparse strata/, test: /^P4 live:/},
  {category: /^Demo studies: live seed$/},
];

export const SHARED_STAND_SKIP_REASON =
  'Writes to the live NCA store; runs only against a local / CI server (localhost, 127.0.0.1)';

/** True for a loopback server — the local stand and the CI stand. */
export function isLoopbackServer(apiRoot: string): boolean {
  let host: string;
  try {
    host = new URL(apiRoot).hostname;
  } catch {
    return false;
  }
  return host === 'localhost' || host === '127.0.0.1' || host === '[::1]';
}

/** True when the rules mark this test as store-writing. */
export function isStoreWriting(category: string, name: string): boolean {
  return STORE_WRITING_TESTS.some((r) => r.category.test(category) && (r.test === undefined || r.test.test(name)));
}

/** The {@link NcaRepository} members that write (rows, snapshots, grants). Reads stay allowed. */
export const STORE_WRITE_METHODS: ReadonlySet<string> = new Set([
  'createStudy', 'renameStudy', 'deleteStudy', 'share', 'createStandalone', 'appendAnalysis',
  'saveAnalysis', 'writeSummary', 'renameAnalysis', 'deleteAnalysis', 'moveToStudy',
  'savePipelinePlan', 'writeSnapshot', 'deleteSnapshot',
]);

export const SHARED_STAND_WRITE_REFUSED =
  'Store write refused on a shared server: add this test to STORE_WRITING_TESTS (tests/shared-stand-guard.ts)';

/**
 * The fail-closed backstop behind the list: on a shared server the tests' repository refuses every
 * write, so a store-writing test missing from {@link STORE_WRITING_TESTS} fails loudly instead of
 * leaving rows behind. Writes that bypass the repository (`ncaDb`, probe users) are covered by the
 * list, and probe users by `withRestrictedUser`'s own refusal.
 */
export function refuseStoreWrites<T extends object>(repo: T): T {
  return new Proxy(repo, {
    get(target, prop, receiver) {
      const v = Reflect.get(target, prop, receiver);
      if (typeof prop === 'string' && STORE_WRITE_METHODS.has(prop) && typeof v === 'function')
        return () => Promise.reject(new Error(`${SHARED_STAND_WRITE_REFUSED} (${prop})`));
      return typeof v === 'function' ? v.bind(target) : v;
    },
  });
}

interface RegisteredTest {readonly name: string; options?: {skipReason?: string}}

/**
 * Mark every store-writing test skipped when the server is not loopback. Runs once, after all test
 * modules registered their categories; the platform's test runner reads `options.skipReason`.
 * Returns how many tests it marked.
 */
export function skipStoreWritingTestsOnSharedStand(
  registry: Record<string, {tests?: RegisteredTest[]}>, apiRoot: string): number {
  if (isLoopbackServer(apiRoot)) return 0;
  let marked = 0;
  for (const [category, entry] of Object.entries(registry)) {
    for (const t of entry.tests ?? []) {
      if (!isStoreWriting(category, t.name)) continue;
      t.options = {...(t.options ?? {}), skipReason: SHARED_STAND_SKIP_REASON};
      marked++;
    }
  }
  return marked;
}
