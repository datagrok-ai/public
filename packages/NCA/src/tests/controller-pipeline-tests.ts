/**
 * Controller — end-to-end pipeline category.
 *
 * @see src/controller/nca-controller.ts (implementation under test)
 *
 * Covers: happy-path completed outcome with results DataFrame and audit
 * entry; validation-failed outcome without compute; engine-throw isolation
 * (one bad subject doesn't blow up the rest); CL label decision; rule-hash
 * stability across object-literal field ordering.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {runNcaPipeline} from '../controller/nca-controller';
import {DEFAULT_RULES, DEFAULT_BLQ_STRATEGY} from '../adapter/defaults';
import {AuditLog, setAuditLog} from '../adapter/audit-log';
import {singleSubjectContract, singleSubjectScalars} from './test-utils';
import type {ColumnContract} from '../adapter/column-contract';
import type {ProfileScalars} from '../adapter/profile-inputs';

function ivBolusContract(): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: 'string'},
      time: {name: 'Time', required: true, rawType: 'float64'},
      concentration: {name: 'conc', required: true, rawType: 'float64'},
    },
  };
}

function ivBolusScalars(): ProfileScalars {
  return {
    dose: 100, route: nca.ROUTE_IV_BOLUS,
    doseUnits: 'mg', concentrationUnits: 'mg/L', timeUnits: 'h',
    bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
  };
}

function cleanIvDf(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject',
      ['S1', 'S1', 'S1', 'S1', 'S1']),
    DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 2, 4, 8])),
    DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 2.5, 0.625, 0.039])),
  ]);
}

category('Controller: end-to-end pipeline', () => {
  test('happy path: completed outcome with 1 result row and 1 audit entry', async () => {
    const log = new AuditLog();
    const outcome = runNcaPipeline(
      cleanIvDf(), ivBolusContract(), ivBolusScalars(), DEFAULT_RULES,
      {log, now: () => 1700000000000},
    );

    expect(outcome.status, 'completed', `status: ${outcome.status}`);
    if (outcome.status !== 'completed') return;
    expect(outcome.results.rowCount, 1, 'expected 1 results row');
    expect(outcome.results.col('SUBJECT')!.get(0), 'S1');
    expect(outcome.results.col('ROUTE')!.get(0), 'IV-bolus');
    // Engine should produce a finite Cmax for the clean profile.
    expect(Number.isFinite(outcome.results.col('CMAX')!.get(0)), true,
      `Cmax must be finite; got ${outcome.results.col('CMAX')!.get(0)}`);

    const entries = log.entries();
    expect(entries.length, 1, 'expected 1 audit entry per subject');
    expect(entries[0].action, 'compute');
    expect(entries[0].subjectKey, 'S1');
    expect(entries[0].timestamp, 1700000000000,
      'options.now injection should drive the audit timestamp');
    expect(entries[0].ruleHash.length > 0, true,
      'ruleHash must be a non-empty fingerprint');
  });

  // GROK-20960 AC-N.4.4. The sci-comp 0.11.0 bump changed what the engine DOES
  // with the declared BLQ rules, but changed no rule VALUE — so the audit
  // trail's fingerprint must be byte-identical across the upgrade, and the
  // `set-zero` the log has always declared is now what the engine actually did.
  // `computeRuleHash` is module-private, so this pins the two things that
  // determine it: the shipped rule values, and the hash of a compute run made
  // with them. A future edit to either breaks this test deliberately.
  test('AC-N.4.4: shipped rule values (and the audit ruleHash) unchanged by the 0.11.0 bump', async () => {
    // (a) The values themselves, spelled out. DEFAULT_BLQ_STRATEGY is what the
    // audit log declares; changing any of these four is a scientific decision,
    // not a refactor.
    expect(DEFAULT_BLQ_STRATEGY.preFirstMeasurable, 'set-zero', 'preFirstMeasurable');
    expect(DEFAULT_BLQ_STRATEGY.embedded, 'set-zero', 'embedded');
    expect(DEFAULT_BLQ_STRATEGY.afterLast, 'set-zero', 'afterLast');
    expect(DEFAULT_BLQ_STRATEGY.consecutiveAfterLast, 'exclude', 'consecutiveAfterLast');
    expect(DEFAULT_RULES.blq, DEFAULT_BLQ_STRATEGY, 'DEFAULT_RULES carries that exact strategy');

    // (b) The fingerprint a compute under those rules writes. Frozen here from
    // the value the controller produced BEFORE the engine bump (the rules are
    // the only input, so an unchanged rule set must reproduce it exactly).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1', 'S1', 'S1']),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Time', [0, 1, 2, 4, 8]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'conc', [0, 10, 7.4, 4.1, 1.2]),
    ]);
    const log = new AuditLog();
    setAuditLog(log);
    const outcome = runNcaPipeline(df, singleSubjectContract(df), singleSubjectScalars(),
      DEFAULT_RULES, {now: () => 1700000000000, log});
    expect(outcome.status, 'completed', 'pipeline must complete');
    const hash = log.entries()[0].ruleHash;
    expect(hash, JSON.stringify({
      aucMethod: DEFAULT_RULES.aucMethod,
      compensatedSummation: DEFAULT_RULES.compensatedSummation,
      extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
      extrapErrorPct: DEFAULT_RULES.extrapErrorPct,
      extrapWarnPctAumc: DEFAULT_RULES.extrapWarnPctAumc,
      blq: {
        preFirstMeasurable: 'set-zero', embedded: 'set-zero',
        afterLast: 'set-zero', consecutiveAfterLast: 'exclude',
      },
      lambdaZ: {
        mode: DEFAULT_RULES.lambdaZ.mode, minPoints: DEFAULT_RULES.lambdaZ.minPoints,
        minRSquared: DEFAULT_RULES.lambdaZ.minRSquared,
        excludeCmax: DEFAULT_RULES.lambdaZ.excludeCmax,
        adjRSquaredFactor: DEFAULT_RULES.lambdaZ.adjRSquaredFactor ?? null,
        // NOTE: `minSpanRatio` is deliberately NOT part of the fingerprint —
        // this package never sets it (the span band is display-side), so it
        // would serialize as a constant. Spelled out because the first cut of
        // this test assumed it was included and the assertion caught that.
        manualPoints: null,
        manualRange: null,
      },
    }), 'the audit ruleHash must be exactly the serialized rule values — no engine-version input');
  });

  test('validation-failed: missing required column → no compute, no audit entries', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
    ]);
    const log = new AuditLog();
    const outcome = runNcaPipeline(
      df, ivBolusContract(), ivBolusScalars(), DEFAULT_RULES, {log},
    );
    expect(outcome.status, 'validation-failed', `expected validation-failed, got ${outcome.status}`);
    if (outcome.status === 'validation-failed') {
      expect(outcome.errors.some((e) => e.code === 'MISSING_COLUMN' && e.column === 'conc'), true,
        `expected MISSING_COLUMN on 'conc'; got ${JSON.stringify(outcome.errors)}`);
    }
    expect(log.entries().length, 0, 'no audit entries on validation-failed');
  });

  test('engine isolation: one bad subject does not block the rest', async () => {
    // S1: clean log-linear decay. S2: all conc = 0 (BLQ everywhere → engine
    // returns status='failed' or throws; controller normalizes to a row of
    // NaN values with status='failed' or 'partial' per spec §6.1 row 3-4).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject',
        ['S1', 'S1', 'S1', 'S1', 'S1', 'S2', 'S2', 'S2']),
      DG.Column.fromFloat64Array('Time',
        Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2])),
      DG.Column.fromFloat64Array('conc',
        Float64Array.from([10, 5, 2.5, 0.625, 0.039, 0, 0, 0])),
    ]);
    const log = new AuditLog();
    const outcome = runNcaPipeline(
      df, ivBolusContract(), ivBolusScalars(), DEFAULT_RULES, {log},
    );
    expect(outcome.status, 'completed');
    if (outcome.status !== 'completed') return;
    expect(outcome.results.rowCount, 2, 'both subjects ship rows');

    // S1 must still ship a finite Cmax — its compute is isolated.
    const statusCol = outcome.results.col('STATUS')!;
    const cmaxCol = outcome.results.col('CMAX')!;
    expect(statusCol.get(0), 'ok',
      `S1 should be ok; got ${statusCol.get(0)} (S2 status: ${statusCol.get(1)})`);
    expect(Number.isFinite(cmaxCol.get(0)), true, 'S1 Cmax must be finite');

    // S2 row exists; status differentiates the two; values NaN per §6.1 row 3.
    const s2Status = statusCol.get(1) as string;
    expect(s2Status === 'failed' || s2Status === 'partial', true,
      `S2 status must be failed|partial; got '${s2Status}'`);

    expect(log.entries().length, 2, 'one audit entry per subject (including failed)');
  });

  test('rule-hash stability across object-literal field reordering', async () => {
    // Two NcaRules instances built with fields in different orders should
    // produce the same audit-log ruleHash — the controller's stable
    // serializer pins the order before stringifying.
    const log1 = new AuditLog();
    const log2 = new AuditLog();
    const rulesA: nca.NcaRules = {
      aucMethod: 'linear-up-log-down',
      compensatedSummation: false,
      blq: {
        preFirstMeasurable: 'set-zero', embedded: 'set-zero',
        afterLast: 'set-zero', consecutiveAfterLast: 'exclude',
      },
      lambdaZ: {
        mode: 'auto-best-fit', minPoints: 3, minRSquared: 0.85,
        excludeCmax: true, adjRSquaredFactor: 1e-4,
      },
      extrapWarnPct: 20, extrapErrorPct: 50, extrapWarnPctAumc: 20,
    };
    const rulesB: nca.NcaRules = {
      extrapErrorPct: 50, extrapWarnPct: 20, extrapWarnPctAumc: 20,
      lambdaZ: {
        adjRSquaredFactor: 1e-4, excludeCmax: true,
        minRSquared: 0.85, minPoints: 3, mode: 'auto-best-fit',
      },
      blq: {
        consecutiveAfterLast: 'exclude', afterLast: 'set-zero',
        embedded: 'set-zero', preFirstMeasurable: 'set-zero',
      },
      compensatedSummation: false, aucMethod: 'linear-up-log-down',
    };

    const o1 = runNcaPipeline(cleanIvDf(), ivBolusContract(), ivBolusScalars(),
      rulesA, {log: log1, now: () => 0});
    const o2 = runNcaPipeline(cleanIvDf(), ivBolusContract(), ivBolusScalars(),
      rulesB, {log: log2, now: () => 0});
    expect(o1.status, 'completed');
    expect(o2.status, 'completed');
    const h1 = log1.entries()[0].ruleHash;
    const h2 = log2.entries()[0].ruleHash;
    expect(h1, h2, `ruleHash must be invariant under key reordering; '${h1}' vs '${h2}'`);
  });

  test('audit entry payload carries the ProfileProvenance snapshot', async () => {
    const log = new AuditLog();
    const outcome = runNcaPipeline(
      cleanIvDf(), ivBolusContract(), ivBolusScalars(), DEFAULT_RULES, {log},
    );
    expect(outcome.status, 'completed');
    if (outcome.status !== 'completed') return;
    const entry = log.entries()[0];
    expect(entry.payload !== undefined && entry.payload !== null, true,
      'audit payload must carry provenance');
    const payload = entry.payload as nca.ProfileProvenance;
    expect(typeof payload.aucMethod, 'string', 'payload.aucMethod must be a string');
    expect(Array.isArray(payload.warnings), true, 'payload.warnings must be an array');
    // The engine should fit λz on a clean log-linear decay; smoke-test the
    // provenance carries it through.
    expect(payload.lambdaZ !== null, true,
      `lambdaZ provenance should be populated on clean profile; got ${JSON.stringify(payload.lambdaZ)}`);
  });

  test('per-profile results: each subject ships independent perProfile slot', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject',
        ['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B']),
      DG.Column.fromFloat64Array('Time',
        Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8])),
      DG.Column.fromFloat64Array('conc',
        Float64Array.from([10, 5, 2.5, 0.625, 0.039, 20, 10, 5, 1.25, 0.078])),
    ]);
    const outcome = runNcaPipeline(df, ivBolusContract(), ivBolusScalars(), DEFAULT_RULES);
    expect(outcome.status, 'completed');
    if (outcome.status !== 'completed') return;
    expect(outcome.perProfile.length, 2);
    expectArray(outcome.perProfile.map((p) => p.index.subjectKey), ['A', 'B']);
  });
});
