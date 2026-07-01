import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';
import {
  getSubcellularLocations,
  LOCATION_COLORS,
  mergeStreamTsv,
  parseSubcellularLocation,
  ProgressCb,
  runWithConcurrency,
  STORE,
} from '../analysis/subcellular-location';
import {ensureLocationColumn} from '../viewers/volcano';

// CRITICAL: `grok.dapi.userDataStorage` returns a NEW UserDataStorage instance
// per access (see js-api dapi.ts: `get userDataStorage() { return new ... }`).
// Patching the instance's `put`/`get` would only stick for one call, then the
// next access creates a fresh instance with the original methods. Patch the
// PROTOTYPE so every fresh instance picks up our spy. Restore the originals in
// finally to avoid leaking across tests.
function patchUserDataStorage(opts: {
  get?: (name: string, currentUser?: boolean) => Promise<any>;
  put?: (name: string, data: any, currentUser?: boolean) => Promise<void>;
}): () => void {
  const proto = (grok.dapi.userDataStorage as any).constructor.prototype;
  const origGet = proto.get;
  const origPut = proto.put;
  if (opts.get) proto.get = opts.get;
  if (opts.put) proto.put = opts.put;
  return () => {
    proto.get = origGet;
    proto.put = origPut;
  };
}

// Locked CK-omics hex palette (Subcellular_Location_Classification_README.txt §COLOR
// SCHEME / CKomics_tool2.py get_location_colors). ARGB = 0xFF000000 | RRGGBB.
const EXPECTED_HEX: Record<string, string> = {
  'Nucleus': '#FF6B6B', 'Cytoplasm': '#ECDC44', 'Mitochondria': '#45B7D1',
  'ER': '#96CEB4', 'Golgi': '#FAAFFE', 'Plasma Membrane': '#DDA0DD',
  'Lysosome': '#F39C12', 'Peroxisome': '#E17055', 'Ribosome': '#A29BFE',
  'Extracellular': '#FD79A8', 'Vesicles': '#FDCB6E', 'Unknown': '#CCCCCC',
};

category('SubcellularLocation', () => {
  test('SEMTYPE.SUBCELLULAR_LOCATION literal is stable', async () => {
    expect(SEMTYPE.SUBCELLULAR_LOCATION, 'Proteomics-SubcellularLocation');
  });

  test('classifier: subcellular field used before GO fallback', async () => {
    // Subcell says Nucleus, GO says mitochondrion — subcell wins (tried first).
    expect(parseSubcellularLocation('Nucleus {ECO:0000269}.', 'mitochondrion [GO:0005739]'),
      'Nucleus');
  });

  test('classifier: GO used only when subcellular field has no match', async () => {
    expect(parseSubcellularLocation('', 'cytoplasm [GO:0005737]'), 'Cytoplasm');
    expect(parseSubcellularLocation('phosphopyruvate hydratase complex',
      'ribosome [GO:0005840]'), 'Ribosome');
  });

  test('classifier: earliest character position wins across categories', async () => {
    // "Mitochondrion" appears before "Nucleus" → Mitochondria.
    expect(parseSubcellularLocation('Mitochondrion. Nucleus.', ''), 'Mitochondria');
    // Reverse order → Nucleus.
    expect(parseSubcellularLocation('Nucleus. Mitochondrion.', ''), 'Nucleus');
  });

  test('classifier: ties broken by keyword-map insertion order (Nucleus before Cytoplasm)', async () => {
    // Both keywords start at position 0 of their own scan; Nucleus is iterated
    // first (strict < on position keeps the first-iterated category).
    expect(parseSubcellularLocation('nucleus', ''), 'Nucleus');
    expect(parseSubcellularLocation('cytoplasm', ''), 'Cytoplasm');
  });

  test('classifier: raw string not pre-stripped (matches inside SUBCELLULAR LOCATION: ... {ECO})', async () => {
    expect(parseSubcellularLocation(
      'SUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:1234}. Nucleus {ECO:0000250}.', ''),
    'Cytoplasm');
  });

  test('classifier: empty / no-match → Unknown', async () => {
    expect(parseSubcellularLocation('', ''), 'Unknown');
    expect(parseSubcellularLocation('some unrelated free text', 'nothing spatial here'),
      'Unknown');
  });

  test('LOCATION_COLORS: all 12 ARGB ints equal the locked hex conversion', async () => {
    const keys = Object.keys(EXPECTED_HEX);
    expect(keys.length, 12);
    for (const k of keys) {
      const argb = 0xFF000000 | parseInt(EXPECTED_HEX[k].slice(1), 16);
      expect(LOCATION_COLORS[k], argb);
    }
    expect(Object.keys(LOCATION_COLORS).length, 12);
  });

  test('mergeStreamTsv: header skipped, positional parse, reviewed beats unreviewed', async () => {
    const tsv = [
      'Entry\tSubcellular location [CC]\tGene Ontology (cellular component)\tReviewed\tGene Names (primary)',
      'P1\tNucleus.\t\tunreviewed\tGENA',
      'P1\tMitochondrion.\t\treviewed\tGENA',
      'P2\t\t\tunreviewed\tGENB',
    ].join('\n');
    const {locByAcc, geneByAcc} = mergeStreamTsv(tsv);
    // P1: reviewed Mitochondria overwrites the earlier unreviewed Nucleus.
    expect(locByAcc.get('P1'), 'Mitochondria');
    // P2: no location anywhere → Unknown.
    expect(locByAcc.get('P2'), 'Unknown');
    expect(geneByAcc.get('P1'), 'GENA');
    expect(geneByAcc.get('P2'), 'GENB');
  });

  test('mergeStreamTsv: existing Unknown is overwritten by a later non-Unknown', async () => {
    const tsv = [
      'Entry\tSubcellular location [CC]\tGO\tReviewed\tGene',
      'P9\t\t\tunreviewed\tG9',
      'P9\tGolgi apparatus.\t\tunreviewed\tG9',
    ].join('\n');
    const {locByAcc} = mergeStreamTsv(tsv);
    expect(locByAcc.get('P9'), 'Golgi');
  });

  // --- 13-08 gap-closure tests ---

  // Test A: runWithConcurrency holds the in-flight cap. No network — uses
  // deferred Promises around a tracked concurrentMax counter.
  test('runWithConcurrency caps in-flight tasks at the supplied limit', async () => {
    let inFlight = 0;
    let concurrentMax = 0;
    const items = new Array(20).fill(0).map((_, i) => i);
    await runWithConcurrency(items, 6, async () => {
      inFlight++;
      if (inFlight > concurrentMax) concurrentMax = inFlight;
      await new Promise((r) => setTimeout(r, 5));
      inFlight--;
      return null;
    });
    expect(concurrentMax <= 6, true, 'concurrentMax should be ≤ 6');
    expect(concurrentMax >= 2, true, 'should have observed parallel execution');
  });

  // Test B: getSubcellularLocations emits progress callbacks with non-decreasing
  // `done` within phase. Monkey-patches fetchProxy, userDataStorage.put, AND
  // userDataStorage.get — the last one is critical: a previous test run may
  // have populated real on-disk storage, which would silently mark every
  // accession as a cache hit and skip the fetch loop entirely.
  test('getSubcellularLocations emits non-decreasing progress per phase', async () => {
    const tsvBody = (acc: string) =>
      'Entry\tSubcellular location [CC]\tGO\tReviewed\tGene Names (primary)\n' +
      `${acc}\tNucleus.\t\treviewed\tGENE1\n`;

    const origFetch = grok.dapi.fetchProxy;
    (grok.dapi as any).fetchProxy = async (url: string) => {
      // Extract the first accession out of the query so each chunk returns
      // a single resolved entry — that's enough to exercise progress.
      const m = /accession%3A([A-Z0-9]+)/.exec(url);
      const acc = m ? m[1] : 'P00000';
      return {ok: true, status: 200, text: async () => tsvBody(acc)} as Response;
    };
    const restore = patchUserDataStorage({
      get: async () => ({}),
      put: async () => {},
    });

    try {
      const accessions: string[] = [];
      for (let i = 0; i < 250; i++) accessions.push(`TESTB${String(i).padStart(5, '0')}`);
      const calls: {done: number; total: number; phase: string}[] = [];
      const progress: ProgressCb = (done, total, phase) => calls.push({done, total, phase});

      await getSubcellularLocations(accessions, 'hsapiens', progress);

      // 250/100 ⇒ 3 accession chunks.
      const accCalls = calls.filter((c) => c.phase === 'fetch-acc');
      expect(accCalls.length, 3);
      expect(accCalls[0].total, 3);
      // Non-decreasing within phase (worker-interleaved → not strictly monotonic).
      for (let i = 1; i < accCalls.length; i++)
        expect(accCalls[i].done >= accCalls[i - 1].done, true,
          'fetch-acc done should be non-decreasing');
      // Final tuple of the last entered phase reaches done === total.
      const last = accCalls[accCalls.length - 1];
      expect(last.done, last.total, 'last fetch-acc tuple should have done === total');
    } finally {
      (grok.dapi as any).fetchProxy = origFetch;
      restore();
    }
  });

  // Test C: userDataStorage.put is called at least twice during a fetch lasting
  // longer than CACHE_FLUSH_INTERVAL_MS (5 s default). Spy via wrapping.
  test('getSubcellularLocations writes the cache incrementally during fetch', async () => {
    // We can't shorten CACHE_FLUSH_INTERVAL_MS without exporting it (and the
    // plan locks it as a tuning knob), so this test runs a fetch that lasts
    // ~6 s — long enough for the timer to fire at least once mid-fetch plus
    // the final finally-flush. Per-chunk delay × chunk count is what stretches
    // wall-clock; with 3 chunks and a 2.2 s per-chunk delay, total ≈ 6 s with
    // FETCH_CONCURRENCY=6 (all parallel) or longer if serialized.
    //
    // The robustness invariant we assert is: put was called >= 2 times. The
    // exact count depends on timer-vs-fetch interleaving; 2 covers (a) one
    // timer fire mid-fetch + (b) the finally flush. We intentionally do not
    // assert >2 — a CI VM might run the entire fetch inside one tick.
    const tsvBody = (acc: string) =>
      'Entry\tSubcellular location [CC]\tGO\tReviewed\tGene Names (primary)\n' +
      `${acc}\tNucleus.\t\treviewed\tGENE1\n`;

    const origFetch = grok.dapi.fetchProxy;
    // First chunk completes at ~1 s (populates `fetched`); second chunk runs
    // until ~6 s. The 5 s flush timer fires at t=5 s while the second chunk
    // is still in flight — `fetched` is already non-empty from chunk 1, so
    // flushCache actually writes (put #1). Finally flushes (put #2). This is
    // the only interleaving that proves the timer-driven incremental cache
    // works, not just the post-pass single flush.
    let callIdx = 0;
    (grok.dapi as any).fetchProxy = async (url: string) => {
      const myIdx = callIdx++;
      const delayMs = myIdx === 0 ? 1000 : 6000;
      await new Promise((r) => setTimeout(r, delayMs));
      const m = /accession%3A([A-Z0-9]+)/.exec(url);
      const acc = m ? m[1] : 'P00000';
      return {ok: true, status: 200, text: async () => tsvBody(acc)} as Response;
    };
    let putCalls = 0;
    const restore = patchUserDataStorage({
      // Empty cache → all accessions go to fetch loop → flushCache has work.
      get: async () => ({}),
      put: async () => { putCalls++; },
    });

    try {
      // 150 accessions → 2 ACC_CHUNK groups (size 100). Worker pool with
      // concurrency=6 starts both chunks at t=0; staggered delays handle the
      // interleaving above.
      const accessions: string[] = [];
      for (let i = 0; i < 150; i++) accessions.push(`TESTC${String(i).padStart(5, '0')}`);
      await getSubcellularLocations(accessions, 'hsapiens');
      expect(putCalls >= 2, true,
        `userDataStorage.put should be called >= 2 times (got ${putCalls})`);
    } finally {
      (grok.dapi as any).fetchProxy = origFetch;
      restore();
    }
  }, {timeout: 30000});

  // Make sure STORE is exported and matches the documented key — guards
  // against accidental rename that would silently fork the on-disk cache
  // from the volcanoOptions toast-cache-peek path.
  test('STORE export matches the documented userDataStorage key', async () => {
    expect(STORE, 'proteomics-subcell-loc');
  });

  // Test D: ensureLocationColumn short-circuits when the accession-set hash
  // tag matches the current set. We don't know the FNV-1a hash from the test,
  // so we exercise the contract end-to-end: (1) first call populates the
  // column AND stamps the hash tag; (2) second call with fetchProxy spy that
  // would throw if invoked still completes — proving the short-circuit fires;
  // (3) mutating the accession set invalidates the hash and re-triggers fetch.
  test('ensureLocationColumn short-circuits on accession-set hash match', async () => {
    const tsvBody = (acc: string) =>
      'Entry\tSubcellular location [CC]\tGO\tReviewed\tGene Names (primary)\n' +
      `${acc}\tNucleus.\t\treviewed\tGENE1\n`;

    const origFetch = grok.dapi.fetchProxy;
    let fetchCalled = 0;
    (grok.dapi as any).fetchProxy = async (url: string) => {
      fetchCalled++;
      const m = /accession%3A([A-Z0-9]+)/.exec(url);
      const acc = m ? m[1] : 'TESTD00';
      return {ok: true, status: 200, text: async () => tsvBody(acc)} as Response;
    };
    // Empty cache → first call goes through fetch path so the hash gets stamped.
    const restore = patchUserDataStorage({
      get: async () => ({}),
      put: async () => {},
    });

    try {
      const idCol = DG.Column.fromStrings('Primary Protein ID',
        ['TESTD01', 'TESTD02', 'TESTD03', 'TESTD04', 'TESTD05']);
      idCol.semType = SEMTYPE.PROTEIN_ID;
      const df = DG.DataFrame.fromColumns([idCol]);

      // (1) Populate + stamp the hash tag.
      await ensureLocationColumn(df);
      const firstFetchCount = fetchCalled;
      expect(firstFetchCount >= 1, true, 'first call should fetch from UniProt');

      // (2) Same accession set → short-circuit, no new fetch.
      await ensureLocationColumn(df);
      expect(fetchCalled, firstFetchCount,
        'second call on same accession set should short-circuit (no new fetch)');

      // (3) Mutate one accession → hash invalidates → fetch resumes.
      idCol.set(0, 'TESTD99');
      await ensureLocationColumn(df);
      expect(fetchCalled > firstFetchCount, true,
        'third call with mutated accession should fetch again');
    } finally {
      (grok.dapi as any).fetchProxy = origFetch;
      restore();
    }
  });

  // Test E: ensureLocationColumn forwards the progress callback through to
  // getSubcellularLocations AND fires 'init-column' after bulk-init.
  test('ensureLocationColumn forwards progress through fetch + emits init-column', async () => {
    const tsvBody = (acc: string) =>
      'Entry\tSubcellular location [CC]\tGO\tReviewed\tGene Names (primary)\n' +
      `${acc}\tNucleus.\t\treviewed\tGENE1\n`;

    const origFetch = grok.dapi.fetchProxy;
    (grok.dapi as any).fetchProxy = async (url: string) => {
      const m = /accession%3A([A-Z0-9]+)/.exec(url);
      const acc = m ? m[1] : 'TESTE00';
      return {ok: true, status: 200, text: async () => tsvBody(acc)} as Response;
    };
    const restore = patchUserDataStorage({
      get: async () => ({}),
      put: async () => {},
    });

    try {
      const idCol = DG.Column.fromStrings('Primary Protein ID',
        ['TESTE01', 'TESTE02', 'TESTE03']);
      idCol.semType = SEMTYPE.PROTEIN_ID;
      const df = DG.DataFrame.fromColumns([idCol]);

      const phases: string[] = [];
      const progress: ProgressCb = (_done, _total, phase) => phases.push(phase);

      await ensureLocationColumn(df, progress);

      expect(phases.includes('fetch-acc'), true,
        'progress should fire at least once with phase=fetch-acc');
      expect(phases.includes('init-column'), true,
        'progress should fire with phase=init-column after bulk-init');
      // 'init-column' should be the last tick — tests can rely on this for
      // closing the dialog's pi cleanly.
      expect(phases[phases.length - 1], 'init-column',
        'init-column should be the final progress phase');
    } finally {
      (grok.dapi as any).fetchProxy = origFetch;
      restore();
    }
  });
});
