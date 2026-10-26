// caMPO spec §2.5 and §10 (questions 1, 2, 4): what Correlation plot exposes and when it recomputes.
// Run as a JavaScript script in Datagrok (Browse | Platform | Functions | Scripts | New JavaScript Script) or paste
// into the console. The results open as the "caMPO corr checks" table and are printed to the browser console.
// For question 4, close the "Correlations (close me)" tab with × afterwards and read the console.

(async () => {
  const report = [];
  const add = (check, ok, details) =>
    report.push({check, result: ok == null ? 'INFO' : ok ? 'PASS' : 'FAIL', details: String(details)});
  const frames = (k = 2) => new Promise((res) => {
    const step = () => (k-- <= 0 ? res() : requestAnimationFrame(step));
    step();
  });
  const tryRun = async (check, fn) => {
    try { await fn(); } catch (e) { add(check, false, `exception: ${e?.message ?? e}`); console.error(e); }
  };

  // ---------- data: desirability-like columns in [0, 1] ----------
  function rng(seed) {
    return () => { seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0; return seed / 4294967296; };
  }
  const gauss = (r) => Math.sqrt(-2 * Math.log(Math.max(r(), 1e-12))) * Math.cos(2 * Math.PI * r());
  const sigm = (x) => 1 / (1 + Math.exp(-x));

  /** A; B, C correlated with A; D independent; E anti-correlated with D; K constant. NaN = missing. */
  function makeData(n, seed, nullPct = 0, extra = 0) {
    const r = rng(seed);
    const names = ['A', 'B', 'C', 'D', 'E', 'K'];
    for (let k = 0; k < extra; k++) names.push(`X${k}`);
    const arr = names.map(() => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      const z0 = gauss(r); const z1 = gauss(r);
      arr[0][i] = sigm(z0);
      arr[1][i] = sigm(z0 + 0.2 * gauss(r));
      arr[2][i] = sigm(z0 + 0.6 * gauss(r));
      arr[3][i] = sigm(z1);
      arr[4][i] = sigm(-0.5 * z1 + gauss(r));
      arr[5][i] = 1;
      for (let k = 0; k < extra; k++) arr[6 + k][i] = sigm(0.3 * z0 + gauss(r));
    }
    if (nullPct > 0) {
      for (let k = 0; k < 5; k++)
        for (let i = 0; i < n; i++) if (r() < nullPct) arr[k][i] = NaN;
    }
    return {names, arr, n};
  }

  function toDf(name, data) {
    const cols = data.names.map((nm, k) =>
      DG.Column.fromFloat32Array(nm, Float32Array.from(data.arr[k], (v) => (Number.isNaN(v) ? DG.FLOAT_NULL : v))));
    const df = DG.DataFrame.fromColumns(cols);
    df.name = name;
    return df;
  }

  // ---------- reference Pearson ----------
  function pearson(a, b, use) {
    let n = 0; let sa = 0; let sb = 0;
    for (let i = 0; i < a.length; i++) if (use(i)) { n++; sa += a[i]; sb += b[i]; }
    const ma = sa / n; const mb = sb / n;
    let sab = 0; let saa = 0; let sbb = 0;
    for (let i = 0; i < a.length; i++) {
      if (!use(i)) continue;
      const x = a[i] - ma; const y = b[i] - mb;
      sab += x * y; saa += x * x; sbb += y * y;
    }
    return sab / Math.sqrt(saa * sbb);
  }

  /** mode: 'listwise' (rows complete in every column) | 'pairwise' (rows complete in both columns). */
  function refMatrix(data, names, mode, rowFilter = null) {
    const idx = names.map((nm) => data.names.indexOf(nm));
    const m = names.length;
    const complete = new Uint8Array(data.n).fill(1);
    for (const k of idx) for (let i = 0; i < data.n; i++) if (Number.isNaN(data.arr[k][i])) complete[i] = 0;
    const R = new Float64Array(m * m);
    for (let p = 0; p < m; p++) {
      for (let q = 0; q < m; q++) {
        const a = data.arr[idx[p]]; const b = data.arr[idx[q]];
        const use = mode === 'listwise' ?
          (i) => complete[i] && (!rowFilter || rowFilter(i)) :
          (i) => !Number.isNaN(a[i]) && !Number.isNaN(b[i]) && (!rowFilter || rowFilter(i));
        R[p * m + q] = pearson(a, b, use);
      }
    }
    return R;
  }

  function viewerMatrix(corr, df, names) {
    const m = names.length;
    const R = new Float64Array(m * m);
    for (let p = 0; p < m; p++)
      for (let q = 0; q < m; q++) R[p * m + q] = corr.getCorrelation(df.col(names[p]), df.col(names[q]));
    return R;
  }

  function maxDiff(R1, R2, m, offDiagonalOnly = true) {
    let d = 0;
    for (let p = 0; p < m; p++) {
      for (let q = 0; q < m; q++) {
        if (offDiagonalOnly && p === q) continue;
        d = Math.max(d, Math.abs(R1[p * m + q] - R2[p * m + q]));
      }
    }
    return d;
  }

  function asymmetry(R, m) {
    let d = 0;
    for (let p = 0; p < m; p++) for (let q = 0; q < p; q++) d = Math.max(d, Math.abs(R[p * m + q] - R[q * m + p]));
    return d;
  }

  const fmtM = (R, m, names) => names.map((nm, p) =>
    nm + ': ' + Array.from({length: m}, (_, q) => R[p * m + q]?.toFixed(3)).join(' ')).join('\n');

  const N = 2000;
  const ACTIVE = ['A', 'B', 'C', 'D', 'E'];
  const m = ACTIVE.length;
  const opts = (names) => ({xColumnNames: names, yColumnNames: names, correlationType: 'Pearson'});

  // ---------- 0. FLOAT_NULL round-trip through Float32Array (spec §2.3) ----------
  const dataNulls = makeData(N, 2, 0.05);
  const dfNulls = toDf('nulls', dataNulls);
  await tryRun('§2.3 FLOAT_NULL in Float32Array', async () => {
    const expected = Array.from(dataNulls.arr[0]).filter(Number.isNaN).length;
    const col = dfNulls.col('A');
    const raw = col.getRawData();
    let eq = 0;
    for (let i = 0; i < raw.length; i++) if (raw[i] === DG.FLOAT_NULL) eq++;
    add('§2.3 missing values recognised by the column', col.stats.missingValueCount === expected,
      `expected ${expected}, stats.missingValueCount = ${col.stats.missingValueCount}`);
    add('§2.3 raw[i] === DG.FLOAT_NULL finds the missing values', eq === expected,
      `raw type ${raw.constructor.name}, matched ${eq} of ${expected}`);
  });

  // ---------- 1. What the viewer exposes (question 2) ----------
  const data0 = makeData(N, 1);
  const df0 = toDf('campo-check', data0);
  const detached = DG.Viewer.correlationPlot(df0, opts(ACTIVE));
  await tryRun('Q2 viewer surface', async () => {
    const tableProp = detached.getOptions(true).look?.table;
    add('Q2 look.table', null, `typeof = ${typeof tableProp}, value = ${JSON.stringify(tableProp)}`);
    const anyTable = detached.table ?? detached.dart?.table;
    add('Q2 JS field `table`', null, anyTable === undefined ? 'undefined' : `typeof = ${typeof anyTable}`);
    add('Q2 getCorrelation exists', typeof detached.getCorrelation === 'function', typeof detached.getCorrelation);
    console.log('[campo] corr plot options', detached.getOptions(true));
  });

  // ---------- 2. Detached viewer (never added to a view): values vs reference ----------
  await tryRun('Q1 detached viewer', async () => {
    const R = viewerMatrix(detached, df0, ACTIVE);
    const ref = refMatrix(data0, ACTIVE, 'listwise');
    const diff = maxDiff(R, ref, m);
    add('Q1 detached: getCorrelation == reference Pearson', diff < 1e-6, `max |Δ| = ${diff.toExponential(2)}`);
    add('Q2 symmetric', asymmetry(R, m) < 1e-9, `max asymmetry = ${asymmetry(R, m).toExponential(2)}`);
    add('Q2 diagonal', null, `getCorrelation(A, A) = ${R[0]}`);
    const k = detached.getCorrelation(df0.col('A'), df0.col('K'));
    add('§2.4 constant column K (not in xColumnNames)', null, `getCorrelation(A, K) = ${k}`);
    console.log('[campo] detached R\n' + fmtM(R, m, ACTIVE));
  });

  // ---------- 3. Missing values: listwise or pairwise deletion? ----------
  await tryRun('§2.4 missing values', async () => {
    const corr = DG.Viewer.correlationPlot(dfNulls, opts(ACTIVE));
    const R = viewerMatrix(corr, dfNulls, ACTIVE);
    const dl = maxDiff(R, refMatrix(dataNulls, ACTIVE, 'listwise'), m);
    const dp = maxDiff(R, refMatrix(dataNulls, ACTIVE, 'pairwise'), m);
    add('§2.4 with 5% missing values', null,
      `vs listwise max |Δ| = ${dl.toExponential(2)}, vs pairwise max |Δ| = ${dp.toExponential(2)}`);
  });

  // ---------- 4. Spearman option is honoured ----------
  await tryRun('§3.2 correlationType', async () => {
    const corr = DG.Viewer.correlationPlot(df0, {...opts(ACTIVE), correlationType: 'Spearman'});
    const s = corr.getCorrelation(df0.col('A'), df0.col('C'));
    const p = detached.getCorrelation(df0.col('A'), df0.col('C'));
    add('§3.2 correlationType: Spearman differs from Pearson', Math.abs(s - p) > 1e-6,
      `Spearman ${s.toFixed(4)}, Pearson ${p.toFixed(4)}`);
  });

  // ---------- 5. Viewer in an inactive FILL tab (question 1) ----------
  const view = grok.shell.addTableView(df0);
  const dm = view.dockManager;
  const corr = DG.Viewer.correlationPlot(df0, opts(ACTIVE));
  await tryRun('Q1 docked FILL, hidden tab', async () => {
    const gridNode = dm.findNode(view.grid.root);
    dm.dock(corr, DG.DOCK_TYPE.FILL, gridNode, 'Correlations');
    await frames();

    let hidden = false;
    const sub = dm.onPanelHidden.subscribe((el) => { if (el === corr.root || el.contains(corr.root)) hidden = true; });
    try {
      gridNode.parent.container.setActiveChild(gridNode.container);
    } catch (e) {
      console.warn('[campo] setActiveChild failed', e);
    }
    await frames();
    sub.unsubscribe();
    const offscreen = corr.root.offsetParent === null;
    add('Q1 corr tab is hidden', hidden || offscreen, `onPanelHidden fired: ${hidden}, offsetParent null: ${offscreen}`);

    // 5a. In-place edit of the bound data frame while hidden.
    const a = df0.col('A');
    for (let i = 0; i < N; i += 2) {
      const v = 1 - data0.arr[0][i];
      data0.arr[0][i] = Math.fround(v);
      a.set(i, v, false);
    }
    a.fireValuesChanged();
    const r1 = viewerMatrix(corr, df0, ACTIVE);
    const d1 = maxDiff(r1, refMatrix(data0, ACTIVE, 'listwise'), m);
    add('Q1 hidden: in-place edit, read immediately', d1 < 1e-6, `max |Δ| = ${d1.toExponential(2)}`);
    await frames(5);
    const d1b = maxDiff(viewerMatrix(corr, df0, ACTIVE), refMatrix(data0, ACTIVE, 'listwise'), m);
    add('Q1 hidden: in-place edit, read after 5 frames', d1b < 1e-6, `max |Δ| = ${d1b.toExponential(2)}`);

    // 5b. Data frame replaced while hidden (what stage A does on every recomputation).
    const data2 = makeData(N, 7);
    const df2 = toDf('campo-check-2', data2);
    const names2 = ['A', 'C', 'D', 'E'];
    corr.dataFrame = df2;
    corr.setOptions(opts(names2));
    const d2 = maxDiff(viewerMatrix(corr, df2, names2), refMatrix(data2, names2, 'listwise'), names2.length);
    add('Q1 hidden: dataFrame replaced, read immediately', d2 < 1e-6, `max |Δ| = ${d2.toExponential(2)}`);

    corr.dataFrame = df0;
    corr.setOptions(opts(ACTIVE));
  });

  // ---------- 6. Row filter of the bound data frame ----------
  await tryRun('row filter', async () => {
    const keep = (i) => i % 3 !== 0;
    df0.filter.init((i) => keep(i));
    await frames();
    const d = maxDiff(viewerMatrix(corr, df0, ACTIVE), refMatrix(data0, ACTIVE, 'listwise', keep), m);
    add('filter: getCorrelation follows df.filter', null,
      `max |Δ| vs filtered reference = ${d.toExponential(2)} (≈0 means the filter is applied)`);
    df0.filter.setAll(true);
  });

  // ---------- 7. Performance (spec §9: n = 1e5, d = 20) ----------
  await tryRun('performance', async () => {
    const big = makeData(100000, 3, 0, 15);
    const names = big.names.filter((nm) => nm !== 'K');
    const dfBig = toDf('big', big);
    const v = DG.Viewer.correlationPlot(dfBig, opts(names));
    let t = performance.now();
    viewerMatrix(v, dfBig, names);
    const first = performance.now() - t;
    t = performance.now();
    viewerMatrix(v, dfBig, names);
    const second = performance.now() - t;
    add('perf: 1e5 rows × 20 columns, full matrix', null,
      `first read ${first.toFixed(0)} ms, second read ${second.toFixed(0)} ms`);
  });

  // ---------- 8. Close events (question 4): close the tab with × by hand ----------
  await tryRun('Q4 close events', async () => {
    const closeMe = DG.Viewer.correlationPlot(df0, opts(ACTIVE));
    dm.dock(closeMe, DG.DOCK_TYPE.FILL, dm.findNode(view.grid.root), 'Correlations (close me)');
    const s1 = grok.events.onViewerClosed.subscribe((args) => {
      const v = args.args?.viewer;
      console.log('[campo] grok.events.onViewerClosed', args, 'viewer === closeMe:', v?.dart === closeMe.dart,
        'type:', v?.type);
    });
    const s2 = dm.onClosed.subscribe((el) => {
      console.log('[campo] dockManager.onClosed', el, 'el === closeMe.root:', el === closeMe.root);
      setTimeout(() => {
        try {
          const r = closeMe.getCorrelation(df0.col('A'), df0.col('B'));
          console.log('[campo] getCorrelation after close:', r);
        } catch (e) {
          console.log('[campo] getCorrelation after close throws:', e?.message ?? e);
        }
        s1.unsubscribe();
        s2.unsubscribe();
      }, 500);
    });
    add('Q4 close events', null, 'close the "Correlations (close me)" tab with × and read the [campo] console lines');
  });

  // A viewer that was never docked throws on close() ("reading 'rN'"), so the detached ones are just dropped.
  console.table(report);
  grok.shell.addTableView(DG.DataFrame.fromObjects(report)).name = 'caMPO corr checks';
})();
