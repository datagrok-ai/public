// Worker-side DG shim tests.
// All run on the main thread, comparing the shim's LiteDataFrame against a
// real DG.DataFrame built from identical input data.

import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {createWorkerDG, arrowIpcToLite, LiteColumn, LiteDataFrame, compileBody,
  clearCompileCache, _getCompileStats, _setCompileCacheCap} from
  '@datagrok-libraries/compute-utils/webworkers';
import {toFeather} from '@datagrok-libraries/arrow';
import {getErrors} from
  '@datagrok-libraries/compute-utils/function-views/src/fitting/fitting-utils';
import {buildSetup, buildRunSeed, LOSS, WorkerPool} from './imports';
import type {OutputTargetItem, ValueBoundsData} from './imports';
import {rangeBound, formulaBound} from './utils';
import {makeWedgedFunc} from './script-fixtures';
import {_package} from '../../../package-test';

function makeBuildSetupArgs(
  inputBounds: Record<string, ValueBoundsData>,
  variedInputNames: string[],
  fixedInputs: Record<string, any>,
) {
  const outputTargets: OutputTargetItem[] = [{
    propName: 'y',
    type: DG.TYPE.FLOAT,
    target: 0,
  }];
  return {
    sessionId: 1,
    fnSource: 'y = a;',
    paramList: ['a', ...Object.keys(fixedInputs)],
    outputParamNames: ['y'],
    lossType: LOSS.RMSE,
    fixedInputs,
    variedInputNames,
    bounds: inputBounds,
    outputTargets,
    nmSettings: new Map<string, number>(),
  };
}

const shim = createWorkerDG();

function expectColumnDataEqual(lite: LiteColumn, dg: DG.Column, msg: string): void {
  expect(lite.length, dg.length, `${msg}: length differs`);
  expect(lite.type, dg.type as any, `${msg}: type differs`);
  const liteRaw = lite.getRawData() as ArrayLike<unknown>;
  const dgRaw = dg.getRawData() as ArrayLike<unknown>;
  expect(liteRaw.length, dgRaw.length, `${msg}: raw length differs`);
  for (let i = 0; i < liteRaw.length; ++i) {
    const a = liteRaw[i];
    const b = dgRaw[i];
    if (typeof a === 'number' && typeof b === 'number')
      expectFloat(a, b, 1e-12, `${msg}[${i}]`);
    else
      expect(String(a), String(b), `${msg}[${i}]`);
  }
}

category('ComputeUtils: Fitting / Worker DG shim', () => {
  test('shim_lite_vs_dg_int_column', async () => {
    const arr = new Int32Array([1, 2, 3, -7, 0, 999]);
    const lite = shim.Column.fromInt32Array('vals', arr);
    const dg = DG.Column.fromInt32Array('vals', arr);
    expectColumnDataEqual(lite, dg, 'int');
    expectFloat(lite.stats.min, dg.stats.min, 1e-12, 'int min');
    expectFloat(lite.stats.max, dg.stats.max, 1e-12, 'int max');
  });

  test('shim_lite_vs_dg_float_column', async () => {
    const arr = new Float64Array([1.5, -2.25, 3.14159, 0, 1e10]);
    const lite = shim.Column.fromFloat64Array('vals', arr);
    const dg = DG.Column.fromFloat64Array('vals', arr);
    expectColumnDataEqual(lite, dg, 'float');
    expectFloat(lite.stats.min, dg.stats.min, 1e-12, 'float min');
    expectFloat(lite.stats.max, dg.stats.max, 1e-12, 'float max');
  });

  test('shim_lite_vs_dg_string_column', async () => {
    const values = ['alpha', 'beta', 'gamma', '', 'alpha', 'beta'];
    const lite = shim.Column.fromStrings('vals', values);
    const dg = DG.Column.fromStrings('vals', values);
    expect(lite.length, dg.length);
    expect(lite.type, DG.COLUMN_TYPE.STRING);
    for (let i = 0; i < values.length; ++i)
      expect(lite.get(i) as string, dg.get(i) as string, `row ${i}`);
  });

  test('shim_lite_vs_dg_bigint_column', async () => {
    let big = BigInt(1);
    for (let i = 0; i < 40; ++i) big = big * BigInt(2);
    const arr = new BigInt64Array([BigInt(1), big, -big, BigInt(0)]);
    const lite = shim.Column.fromBigInt64Array('vals', arr);
    const dg = DG.Column.fromBigInt64Array('vals', arr);
    expect(lite.length, dg.length, 'length differs');
    expect(lite.type, dg.type as any, 'type differs');
    // DG.Column.fromBigInt64Array.getRawData() throws "Not implemented" on the
    // current DG version; assert against the original input array instead.
    const liteRaw = lite.getRawData() as unknown as BigInt64Array;
    for (let i = 0; i < arr.length; ++i)
      expect(String(liteRaw[i]), String(arr[i]), `row ${i}`);
  });

  test('shim_lite_vs_dg_datetime_column', async () => {
    // Pin parity for datetime columns: get(i) must return a dayjs-like value
    // matching DG's `dayjs(dart)` wrap (js-api/src/wrappers_impl.ts:82-83);
    // getRawData() must be the microsecond Float64Array DG uses; toList()
    // must contain dayjs-like values. Compare via .valueOf() (ms epoch) so
    // the assertion is independent of cross-realm dayjs identity.
    const vals = [
      dayjs('2024-01-01T00:00:00Z'),
      dayjs('2024-06-15T12:30:00Z'),
      dayjs('1970-01-01T00:00:00Z'),
    ];
    const lite = shim.Column.fromList(shim.COLUMN_TYPE.DATE_TIME, 't', vals);
    const dg = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 't', vals);
    expect(lite.length, dg.length, 'length differs');
    expect(lite.type, dg.type as any, 'type differs');
    const liteRaw = lite.getRawData() as Float64Array;
    const dgRaw = dg.getRawData() as Float64Array;
    expect(liteRaw.length, dgRaw.length, 'raw length differs');
    for (let i = 0; i < liteRaw.length; ++i)
      expectFloat(liteRaw[i], dgRaw[i], 1e-12, `raw[${i}]`);
    for (let i = 0; i < lite.length; ++i) {
      const a = lite.get(i) as any; const b = dg.get(i) as any;
      expectFloat(Number(a?.valueOf?.()), Number(b?.valueOf?.()), 1e-9,
        `get(${i}) ms epoch`);
      expect(typeof a?.format, 'function', `lite get(${i}) lacks dayjs.format`);
      expect(typeof b?.format, 'function', `dg get(${i}) lacks dayjs.format`);
    }
    // Lite's toList must be self-consistent with its get(i). DG's
    // Column.toList() for datetime returns nulls (separate DG quirk that
    // doesn't agree with DG.Column.get); we don't assert against it.
    const liteList = lite.toList();
    expect(liteList.length, lite.length, 'toList length');
    for (let i = 0; i < liteList.length; ++i)
      expectFloat(Number((liteList[i] as any)?.valueOf?.()),
        Number((lite.get(i) as any)?.valueOf?.()), 1e-9,
        `toList[${i}] vs get(${i})`);
  });

  test('shim_dataframe_fromcolumns_parity', async () => {
    const ids = new Int32Array([10, 20, 30]);
    const vals = new Float64Array([1.1, 2.2, 3.3]);
    const liteDf = shim.DataFrame.fromColumns([
      shim.Column.fromInt32Array('id', ids),
      shim.Column.fromFloat64Array('value', vals),
      shim.Column.fromStrings('label', ['x', 'y', 'z']),
    ]);
    const dgDf = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('id', ids),
      DG.Column.fromFloat64Array('value', vals),
      DG.Column.fromStrings('label', ['x', 'y', 'z']),
    ]);
    expect(liteDf.rowCount, dgDf.rowCount, 'row count');
    expect(liteDf.columns.length, dgDf.columns.length, 'col count');
    expect(liteDf.columns.names().join(','), dgDf.columns.names().join(','), 'col order');
    expectColumnDataEqual(liteDf.col('id')!, dgDf.col('id')!, 'id col');
    expectColumnDataEqual(liteDf.col('value')!, dgDf.col('value')!, 'value col');
    // string col compared via .get(i)
    const a = liteDf.col('label')!;
    const b = dgDf.col('label')!;
    for (let i = 0; i < a.length; ++i)
      expect(a.get(i) as string, b.get(i) as string);
  });

  test('shim_geterrors_parity', async () => {
    // Build a target DF and a simulation DF, twice — once as DG, once as
    // shim Lite — and run getErrors() on both. Output must be byte-identical.
    const argRaw = new Float64Array([0, 0.5, 1.0, 1.5, 2.0]);
    const fnRaw = new Float64Array([10, 6.07, 3.68, 2.23, 1.35]);
    const simArgRaw = new Float64Array([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]);
    const simFnRaw = new Float64Array([10.1, 7.8, 6.0, 4.7, 3.6, 2.8, 2.2, 1.7, 1.3]);

    // Real DG side
    const expArgDg = DG.Column.fromFloat64Array('t', argRaw);
    const expFnDg = DG.Column.fromFloat64Array('y', fnRaw);
    const simDg = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('t', simArgRaw),
      DG.Column.fromFloat64Array('y', simFnRaw),
    ]);

    // Shim Lite side
    const expArgLt = shim.Column.fromFloat64Array('t', argRaw);
    const expFnLt = shim.Column.fromFloat64Array('y', fnRaw);
    const simLt = shim.DataFrame.fromColumns([
      shim.Column.fromFloat64Array('t', simArgRaw),
      shim.Column.fromFloat64Array('y', simFnRaw),
    ]);

    const dgErrs = getErrors(expArgDg, [expFnDg], simDg, true);
    const ltErrs = getErrors(expArgLt as any, [expFnLt as any], simLt as any, true);
    expect(dgErrs.length, ltErrs.length, 'errors length');
    for (let i = 0; i < dgErrs.length; ++i)
      expectFloat(dgErrs[i], ltErrs[i], 1e-12, `err[${i}]`);
  });

  test('arrow_to_lite_roundtrip_int', async () => {
    const dgDf = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('vals', new Int32Array([1, 2, 3, -7, 0])),
    ]);
    const bytes = toFeather(dgDf, true)!;
    const lite = arrowIpcToLite(bytes);
    expect(lite.rowCount, dgDf.rowCount);
    expectColumnDataEqual(lite.col('vals')!, dgDf.col('vals')!, 'int roundtrip');
  });

  test('arrow_to_lite_roundtrip_float', async () => {
    const dgDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('vals', new Float64Array([1.5, -2.25, 3.14159, 0, 1e10])),
    ]);
    const bytes = toFeather(dgDf, true)!;
    const lite = arrowIpcToLite(bytes);
    expect(lite.rowCount, dgDf.rowCount);
    expectColumnDataEqual(lite.col('vals')!, dgDf.col('vals')!, 'float roundtrip');
  });

  test('arrow_to_lite_roundtrip_string', async () => {
    const values = ['alpha', 'beta', 'gamma', '', 'alpha'];
    const dgDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('vals', values)]);
    const bytes = toFeather(dgDf, true)!;
    const lite = arrowIpcToLite(bytes);
    const c = lite.col('vals')!;
    expect(c.length, values.length);
    for (let i = 0; i < values.length; ++i)
      expect(c.get(i) as string, values[i], `row ${i}`);
  });

  test('arrow_to_lite_roundtrip_datetime', async () => {
    // Catches the arrow-to-lite double-×1000 regression: if the Arrow path
    // doesn't store microseconds, get(i) values come back ~10⁶× off.
    const vals = [
      dayjs('2024-01-01T00:00:00Z'),
      dayjs('2024-06-15T12:30:00Z'),
      dayjs('1970-01-01T00:00:00Z'),
    ];
    const dgDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 't', vals),
    ]);
    const bytes = toFeather(dgDf, true)!;
    const lite = arrowIpcToLite(bytes);
    const liteCol = lite.col('t')!;
    const dgCol = dgDf.col('t')!;
    expect(liteCol.length, dgCol.length, 'length differs');
    expect(liteCol.type, dgCol.type as any, 'type differs');
    const liteRaw = liteCol.getRawData() as Float64Array;
    const dgRaw = dgCol.getRawData() as Float64Array;
    for (let i = 0; i < liteRaw.length; ++i)
      expectFloat(liteRaw[i], dgRaw[i], 1e-12, `raw[${i}]`);
    for (let i = 0; i < liteCol.length; ++i) {
      const a = liteCol.get(i) as any; const b = dgCol.get(i) as any;
      expectFloat(Number(a?.valueOf?.()), Number(b?.valueOf?.()), 1e-9,
        `roundtrip get(${i}) ms epoch`);
      expect(typeof a?.format, 'function', `lite get(${i}) lacks dayjs.format`);
    }
  });

  test('shim_stats_lazy_caching', async () => {
    const arr = new Float64Array([1, 2, 3, 4, 5]);
    const lite = shim.Column.fromFloat64Array('x', arr);
    const s1 = lite.stats;
    const s2 = lite.stats;
    expect(s1 === s2, true, 'second access must hit cache (same object identity)');
    expectFloat(s1.min, 1, 1e-12);
    expectFloat(s1.max, 5, 1e-12);
  });

  test('shim_body_wrap_captures_outputs', async () => {
    // Mimics the worker entry's body-wrap: given a script body that assigns
    // to a declared output variable, compile it with the shim's DG injected
    // and assert the captured output is a LiteDataFrame.
    const body = "simulation = DG.DataFrame.fromColumns([" +
      "DG.Column.fromFloat64Array('y', new Float64Array([1,2,3]))" +
      "]);";
    const wrapped = `var simulation; ${body}; return {simulation};`;
    const fn = new Function('DG', wrapped) as (dg: any) => {simulation: LiteDataFrame};
    const out = fn(shim);
    expect(out.simulation != null, true, 'simulation output captured');
    expect(out.simulation.rowCount, 3);
    const c = out.simulation.col('y')!;
    expect(c.length, 3);
    expect(c.type, 'double');
    const raw = c.getRawData() as Float64Array;
    expectFloat(raw[0], 1, 1e-12);
    expectFloat(raw[2], 3, 1e-12);
  });

  test('shim_realscript_lorenz_smoke', async () => {
    // Load the real production script body, parse the //input: header,
    // compile via the production compileBody with the shim, and assert the
    // output is a LiteDataFrame with the expected columns. This is the
    // sentinel for shim coverage gaps — if a real script uses a DG.* surface
    // the shim doesn't expose, this test fails immediately.
    const src = await _package.files.readAsText('shim-fixtures/lorenz-attractor.js');
    const {body, inputNames, outputNames} = parseScript(src);
    const fn = compileBody(body, inputNames, outputNames);
    const inputs = {iterations: 100, dt: 0.01, x0: 0, y0: 1, z0: 1.05};
    const out = fn(shim, ...inputNames.map((n) => (inputs as any)[n])) as Record<string, LiteDataFrame>;
    expect(out.df != null, true, 'df output captured');
    expect(out.df.rowCount, 101, 'iterations + 1 rows');
    expect(out.df.columns.names().join(','), 'X,Y,Z');
  });

  test('shim_realscript_object_cooling_smoke', async () => {
    // Regression: object-cooling.js declares `const tempDiff` and
    // `const coolingFactor` whose names also appear in //output headers.
    // Predeclaring outputs as `let X` would SyntaxError; the wrap relies on
    // JS scope resolution at the return site instead.
    const src = await _package.files.readAsText('shim-fixtures/object-cooling.js');
    const {body, inputNames, outputNames} = parseScript(src);
    const fn = compileBody(body, inputNames, outputNames);
    // simTime large enough for the cube to actually cool past desiredTemp
    // (~19280s for these defaults), so `timeToCool` gets assigned.
    const inputs = {ambTemp: 22, initTemp: 100, desiredTemp: 30,
      area: 0.06, heatCap: 4200, heatTransferCoeff: 8.3,
      previousRun: null, simTime: 25000};
    const out = fn(shim, ...inputNames.map((n) => (inputs as any)[n])) as Record<string, any>;
    expect(out.simulation != null, true, 'simulation output captured');
    expect(out.simulation.rowCount, 25000);
    expect(out.simulation.columns.names().join(','), 'Time,Temperature');
    expectFloat(out.tempDiff, 78, 1e-12, 'tempDiff = initTemp - ambTemp');
    expectFloat(out.coolingFactor, 8.3 * 0.06 / 4200, 1e-12, 'coolingFactor');
    expect(typeof out.timeToCool, 'number', 'timeToCool was assigned');
  });

  test('shim_realscript_heat_exchange_smoke', async () => {
    const src = await _package.files.readAsText('shim-fixtures/heat-exchange.js');
    const {body, inputNames, outputNames} = parseScript(src);
    const fn = compileBody(body, inputNames, outputNames);
    const inputs = {len: 100, k: 1.2};
    const out = fn(shim, ...inputNames.map((n) => (inputs as any)[n])) as Record<string, any>;
    expect(out.simulation != null, true, 'simulation output captured');
    expect(out.simulation.rowCount, 100);
    expect(out.simulation.columns.names().join(','), 'time,concentration,temp,saturation');
    expect(typeof out.finalTemperature, 'number');
    expect(typeof out.finalSaturation, 'number');
    expect(typeof out.finalConcentration, 'number');
  });

  test('shim_realscript_lotka_volterra_smoke', async () => {
    const src = await _package.files.readAsText('shim-fixtures/lotka-volterra.js');
    const {body, inputNames, outputNames} = parseScript(src);
    const fn = compileBody(body, inputNames, outputNames);
    const inputs = {x0: 0.5, y0: 2, alpha: 0.9, beta: 0.8,
      gamma: 0.9, sigma: 0.7, h: 0.01, timeStart: 0, timeStop: 18};
    const out = fn(shim, ...inputNames.map((n) => (inputs as any)[n])) as Record<string, LiteDataFrame>;
    expect(out.df != null, true, 'df output captured');
    expect(out.df.columns.names().join(','), 'Time,Number of preys,Number of predators');
    expect(out.df.rowCount > 0, true, 'integration produced rows');
  });

  test('pool_dispose_drains_pending_setups', async () => {
    // White-box: dispose's job is to resolve any pendingSetups it finds.
    // Inject one directly so the assertion doesn't depend on Worker
    // round-trip timing (which would only happen to be deterministic
    // because of JS's single-threaded event loop, but is brittle to read).
    type SlotInternal = {
      pendingSetups: Map<number, {sessionId: number; resolve: (ack: SetupAckLike) => void}>;
    };
    type SetupAckLike = {kind: 'setup-ack'; sessionId: number; ok: boolean; message?: string};
    const pool = new WorkerPool(1);
    const internal = pool as unknown as {
      ensureWorkers: () => void;
      slots: SlotInternal[];
    };
    internal.ensureWorkers();
    let resolved = false;
    let ack: SetupAckLike | null = null;
    internal.slots[0].pendingSetups.set(7, {
      sessionId: 7,
      resolve: (a) => { resolved = true; ack = a; },
    });
    pool.dispose();
    expect(resolved, true, 'dispose must call pendingSetups[*].resolve');
    expect(ack !== null, true);
    expect(ack!.kind, 'setup-ack');
    expect(ack!.ok, false, 'drained ack must signal failure');
    expect(/disposed/i.test(ack!.message ?? ''), true,
      `expected disposed-flavor message, got: ${ack!.message}`);
  });

  test('pool_dispose_drains_inflight_run', async () => {
    // White-box: dispose's job is to resolve the in-flight run as a failure.
    // Inject a fake `running` slot state so the test doesn't hinge on the
    // dispatchRun → worker → reply round-trip.
    type RunReplyLike = {kind: string; taskId: number; message?: string};
    type SlotInternal = {
      runState:
        | {phase: 'idle'}
        | {
            phase: 'running';
            run: {
              run: {taskId: number; sessionId: number; kind: string; seed: Float64Array};
              transferables: Transferable[];
              resolve: (r: RunReplyLike) => void;
            };
            runTimer: ReturnType<typeof setTimeout>;
          };
    };
    const pool = new WorkerPool(1);
    const internal = pool as unknown as {
      ensureWorkers: () => void;
      slots: SlotInternal[];
    };
    internal.ensureWorkers();
    let resolved = false;
    let reply: RunReplyLike | null = null;
    // Long-lived noop timer — dispose() will clearTimeout it as part of
    // the running→idle transition.
    const noopTimer = setTimeout(() => {}, 60_000);
    internal.slots[0].runState = {
      phase: 'running',
      run: {
        run: {taskId: 42, sessionId: 9, kind: 'run-seed', seed: new Float64Array([1, 2])},
        transferables: [],
        resolve: (r) => { resolved = true; reply = r; },
      },
      runTimer: noopTimer,
    };
    pool.dispose();
    expect(resolved, true, 'dispose must resolve the running slot\'s run');
    expect(reply !== null, true);
    expect(reply!.kind, 'failure');
    expect(reply!.taskId, 42);
    expect(/disposed/i.test(reply!.message ?? ''), true,
      `expected disposed-flavor message, got: ${reply!.message}`);
  });

  test('compile_cache_hits_across_compiles', async () => {
    // Module-level compile cache must hit when the same body+param+output
    // signature recompiles, so a long-lived worker pays parse cost once
    // per source. Drive compileBody directly on the main thread — same
    // module instance, same cache.
    clearCompileCache();
    const bodyA = 'y = a + 1;';
    const bodyB = 'y = a * 2;';
    compileBody(bodyA, ['a'], ['y']);
    let s = _getCompileStats();
    expect(s.misses, 1, 'first compile is a miss');
    expect(s.hits, 0);
    compileBody(bodyA, ['a'], ['y']);
    compileBody(bodyA, ['a'], ['y']);
    s = _getCompileStats();
    expect(s.misses, 1, 'identical body should not recompile');
    expect(s.hits, 2);
    compileBody(bodyB, ['a'], ['y']);
    s = _getCompileStats();
    expect(s.misses, 2, 'different body should compile fresh');
    compileBody(bodyA, ['a'], ['y']);
    s = _getCompileStats();
    expect(s.hits, 3, 'A still cached after B was added');
    expect(s.misses, 2);
  });

  test('compile_cache_evicts_lru', async () => {
    // Cap the cache at 1 so adding B evicts A; the next A re-compile
    // proves LRU eviction is wired and observable.
    _setCompileCacheCap(1);
    const bodyA = 'y = a + 10;';
    const bodyB = 'y = a + 20;';
    compileBody(bodyA, ['a'], ['y']);
    compileBody(bodyB, ['a'], ['y']);  // evicts A
    compileBody(bodyA, ['a'], ['y']);  // miss → recompile
    const s = _getCompileStats();
    expect(s.misses, 3, 'A→B→A with cap=1 must produce 3 misses');
    expect(s.hits, 0);
    // Restore cap to default for downstream tests.
    _setCompileCacheCap(32);
  });

  test('pool_replacement_reprimes_active_session', async () => {
    // After removeSlot replaces a worker mid-fit, the new slot must be
    // re-primed for any session whose `setupAll` already succeeded so
    // that a subsequent dispatchRun on that session reaches the
    // replacement instead of rejecting with "no slot primed".
    const pool = new WorkerPool(1);
    try {
      const targets: OutputTargetItem[] = [{
        propName: 'y', type: DG.TYPE.FLOAT, target: 0,
      }];
      const inputBounds: Record<string, ValueBoundsData> = {a: rangeBound(0, 5, 'a')};
      const {setup} = buildSetup({
        sessionId: 99, fnSource: 'y = a + 1;',
        paramList: ['a'], outputParamNames: ['y'],
        lossType: LOSS.RMSE, fixedInputs: {}, variedInputNames: ['a'],
        bounds: inputBounds, outputTargets: targets, nmSettings: new Map(),
      });
      await pool.setupAll(setup, []);
      expect(pool._slotsForTest(), 1);

      // Force the only slot out; pool spawns a fresh replacement and
      // re-primes it for session 99.
      const internal = pool as unknown as {
        removeSlot: (slot: any, reason: string) => void;
        slots: any[];
      };
      internal.removeSlot(internal.slots[0], 'test forced remove');
      expect(pool._slotsForTest(), 1, 'replacement must be spawned');

      // dispatchRun on the (now-replaced) session must succeed via the
      // re-primed slot. Pre-fix this would reject synchronously with
      // "no slot is primed" because the replacement has empty sessions.
      const {run, transferables} = buildRunSeed({
        sessionId: 99, taskId: 0, seedIndex: 0, seed: new Float64Array([1]),
      });
      const reply = await pool.dispatchRun({run, transferables});
      expect(reply.kind, 'success', 'replacement must serve the active session');
    } finally {
      pool.dispose();
    }
  });

  test('pool_replaces_slot_on_remove', async () => {
    // White-box: removeSlot must spawn a fresh slot so a long-lived pool
    // doesn't drift toward zero workers under repeated errors/timeouts.
    type SlotInternal = {worker: Worker};
    const pool = new WorkerPool(2);
    const internal = pool as unknown as {
      ensureWorkers: () => void;
      removeSlot: (slot: SlotInternal, reason: string) => void;
      slots: SlotInternal[];
    };
    try {
      internal.ensureWorkers();
      expect(pool._slotsForTest(), 2);
      const originalWorker0 = internal.slots[0].worker;
      internal.removeSlot(internal.slots[0], 'test forced remove');
      expect(pool._slotsForTest(), 2, 'removeSlot must spawn a replacement');
      // The replacement is appended at the end; the surviving original
      // (was at index 1) is now at index 0. Either way, the freshly
      // spawned worker must not be the one we terminated.
      for (const slot of internal.slots)
        expect(slot.worker !== originalWorker0, true,
          'replacement must be a fresh Worker');
    } finally {
      pool.dispose();
    }
  });

  test('pool_run_timeout_unblocks_wedged_worker', async () => {
    // Worker body is `while (true)` — the cost-function call never returns.
    // Without runTimeoutMs the dispatchRun promise hangs forever; with it,
    // the slot is terminated and the run resolves as a timed-out failure.
    const pool = new WorkerPool(1, {runTimeoutMs: 500});
    try {
      const func = makeWedgedFunc();
      const targets: OutputTargetItem[] = [{
        propName: 'y', type: DG.TYPE.FLOAT, target: 0,
      }];
      const inputBounds: Record<string, ValueBoundsData> = {a: rangeBound(0, 5, 'a')};
      const {setup} = buildSetup({
        sessionId: 2, fnSource: (func as DG.Script).script,
        paramList: ['a'], outputParamNames: ['y'],
        lossType: LOSS.RMSE, fixedInputs: {}, variedInputNames: ['a'],
        bounds: inputBounds, outputTargets: targets, nmSettings: new Map(),
      });
      await pool.setupAll(setup, []);  // setup is fast — only compile, no run
      const {run, transferables} = buildRunSeed({
        sessionId: 2, taskId: 0, seedIndex: 0, seed: new Float64Array([1]),
      });
      const start = Date.now();
      const reply = await pool.dispatchRun({run, transferables});
      const elapsed = Date.now() - start;
      expect(reply.kind, 'failure', 'wedged run must resolve as failure');
      expect(/timed? ?out/i.test((reply as any).message), true,
        `expected timeout-flavor message, got: ${(reply as any).message}`);
      expect(elapsed < 2000, true, `expected ~500ms, got ${elapsed}ms`);
    } finally {
      pool.dispose();
    }
  });

  test('bounds_const_dataframe_strips_for_postmessage', async () => {
    // A const DG.DataFrame in inputBounds.value would crash structuredClone
    // before any worker code runs. After buildSetup the value is null;
    // the DataFrame travels via setup.fixedDataFrames as Arrow IPC bytes.
    const refDf = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.INT, 'k', [1, 2, 3])]);
    const inputBounds: Record<string, ValueBoundsData> = {
      refDf: {type: 'const', value: refDf},
      a: rangeBound(0, 5, 'a'),
    };
    const fixedInputs = {refDf};
    const {setup} = buildSetup(makeBuildSetupArgs(inputBounds, ['a'], fixedInputs));
    expect((setup.bounds.refDf as any).type, 'const');
    expect((setup.bounds.refDf as any).value, null);
    expect(setup.fixedDataFrames.refDf instanceof Uint8Array, true, 'DataFrame in fixedDataFrames');
    // The acid test: structured-cloning the setup must not throw.
    structuredClone(setup);
  });

  test('bounds_const_dayjs_strips_for_postmessage', async () => {
    const t0 = dayjs('2024-01-01T00:00:00Z');
    const inputBounds: Record<string, ValueBoundsData> = {
      t0: {type: 'const', value: t0},
      a: rangeBound(0, 5, 'a'),
    };
    const fixedInputs = {t0};
    const {setup} = buildSetup(makeBuildSetupArgs(inputBounds, ['a'], fixedInputs));
    expect((setup.bounds.t0 as any).type, 'const');
    expect((setup.bounds.t0 as any).value, null);
    expect(setup.fixedInputTypes.t0, 'dayjs');
    expect(setup.fixedInputs.t0 as any, t0.toISOString());
    structuredClone(setup);
  });

  test('bounds_const_date_strips_for_postmessage', async () => {
    const t0 = new Date('2024-01-01T00:00:00Z');
    const inputBounds: Record<string, ValueBoundsData> = {
      t0: {type: 'const', value: t0},
      a: rangeBound(0, 5, 'a'),
    };
    const fixedInputs = {t0};
    const {setup} = buildSetup(makeBuildSetupArgs(inputBounds, ['a'], fixedInputs));
    expect((setup.bounds.t0 as any).type, 'const');
    expect((setup.bounds.t0 as any).value, null);
    expect(setup.fixedInputTypes.t0, 'date');
    expect(setup.fixedInputs.t0 as any, t0.toISOString());
    structuredClone(setup);
  });

  test('bounds_changing_entries_passthrough_unchanged', async () => {
    // Sanity: `changing` bounds (numeric or formula) round-trip through
    // buildSetup with all fields intact. Only `const` entries get stripped.
    const inputBounds: Record<string, ValueBoundsData> = {
      a: {type: 'const', value: 5},
      x: formulaBound('a-1', 'a+1', 'x'),
      y: rangeBound(-10, 10, 'y'),
    };
    const {setup} = buildSetup(makeBuildSetupArgs(inputBounds, ['x', 'y'], {a: 5}));
    expect((setup.bounds.a as any).value, null, 'const value stripped');
    const xb = setup.bounds.x as any;
    expect(xb.type, 'changing', 'formula bound type preserved');
    expect(xb.bottom.type, 'formula');
    expect(xb.bottom.formula, 'a-1');
    expect(xb.top.type, 'formula');
    expect(xb.top.formula, 'a+1');
    const yb = setup.bounds.y as any;
    expect(yb.type, 'changing', 'value bound type preserved');
    expect(yb.bottom.type, 'value');
    expect(yb.bottom.value, -10);
    expect(yb.top.type, 'value');
    expect(yb.top.value, 10);
    structuredClone(setup);
  });

});

/** Minimal Datagrok-script header parser: extracts input names and output names. */
function parseScript(src: string): {body: string; inputNames: string[]; outputNames: string[]} {
  const inputNames: string[] = [];
  const outputNames: string[] = [];
  const lines = src.split('\n');
  let bodyStart = 0;
  for (let i = 0; i < lines.length; ++i) {
    const ln = lines[i];
    const trimmed = ln.trim();
    if (trimmed.startsWith('//')) {
      // //input: <type> <name> [= ...] [{...}] | //output: <type> <name> [...]
      const m = trimmed.match(/^\/\/(input|output):\s*\S+\s+([A-Za-z_$][A-Za-z0-9_$]*)/);
      if (m) {
        if (m[1] === 'input') inputNames.push(m[2]);
        else outputNames.push(m[2]);
      }
      bodyStart = i + 1;
    } else if (trimmed === '') {
      bodyStart = i + 1;
    } else {
      break;
    }
  }
  return {
    body: lines.slice(bodyStart).join('\n'),
    inputNames, outputNames,
  };
}
