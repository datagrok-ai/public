// Grokky latency/accuracy benchmark harness — see docs/BENCHMARK.md.
//
// Drives the golden prompt suite (files/benchmark/suite.yaml) through the SAME runtime path the
// panel uses (workspace context + view tools + real datagrok-exec/verify round-trips), times each
// turn client-side (TTFT / total / tool round-trips), reads the server-side token/cost metrics the
// runtime forwards on `final`, scores each prompt (deterministic assert + Haiku judge), and downloads
// a JSON + Markdown report tagged with a config label so runs are comparable across latency levers.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as jsyaml from 'js-yaml';

import {_package} from '../../package';
import {ClaudeRuntimeClient, ClaudeModel, TurnServerMetrics} from '../../claude/runtime-client';
import {viewFunctionTools} from '../view-tools';
import {executeSingleBlock, runVerification, buildWorkspaceContext} from '../../claude/exec-blocks';

interface SuitePrompt {
  category: string;
  prompt: string;
  table?: string;
  assert?: string;
  rubric?: string;
  /** Drives the per-difficulty rollup, so a comparison can be read per tier. */
  difficulty?: 'trivial' | 'standard' | 'hard';
}

interface RepTurn {
  ok: boolean;
  ttftMs: number | null;
  totalMs: number | null;
  toolRoundTrips: number;
  /** Bare tool names in call order — lets asserts check WHICH path a turn took. */
  tools: string[];
  ctxChars: number;
  toolDefCount: number;
  server: TurnServerMetrics | null;
  assertPass: boolean | null;
  /** The runtime's own "could not confirm this action took effect" flag. */
  unverified: boolean;
  content: string;
  error?: string;
  /** Set when an assert expression itself threw — distinguishes a broken assert from a model failure. */
  assertError?: string;
  /** Compact description of the workspace when an assert failed, so the report says what the model
   * actually did rather than only that it was wrong. */
  failState?: string;
}

interface Stats {n: number; median: number | null; mean: number | null; min: number | null; max: number | null; p90: number | null; std: number | null}
interface JudgeResult {pass: boolean; score: number; reason: string}

interface PromptResult {
  category: string;
  prompt: string;
  difficulty: 'trivial' | 'standard' | 'hard';
  /** The model this arm pinned, or 'default' when the runtime chose. */
  model: string;
  reps: number;
  okCount: number;
  /** Union of tool names seen across reps — the cheap-vs-strong path discriminator. */
  toolsUsed: string[];
  unverifiedRate: number;
  ttftMs: Stats;
  totalMs: Stats;
  toolRoundTrips: Stats;
  inputTokens: Stats;
  outputTokens: Stats;
  cacheReadTokens: Stats;
  costUsd: Stats;
  ctxChars: number;
  toolDefCount: number;
  assertPassRate: number | null;
  judge: JudgeResult | null;
  accuracyPass: boolean | null;
  /** Workspace state from the first failing rep — what the model actually left behind. */
  failState?: string;
  errors: string[];
}

interface BenchReport {label: string; when: string; reps: number; model: string; hasTokens: boolean;
  /** The `only` filter this arm ran with, so a slice is never mistaken for a full run. */
  only?: string;
  results: PromptResult[]}

// A hard prompt (substructure search, multi-table chemistry) legitimately runs past two minutes.
// Timing it out scores a slow turn as a *wrong* one, which conflates latency with correctness —
// the two axes this harness exists to separate. Latency is still measured; only the cutoff moves.
const TURN_TIMEOUT_MS: Record<string, number> = {trivial: 120000, standard: 120000, hard: 240000};
const turnTimeout = (difficulty?: string): number => TURN_TIMEOUT_MS[difficulty ?? 'standard'] ?? 120000;
// A UI FilterGroup filter applies asynchronously — `t.filter` updates a frame or two after
// updateOrAdd (see the datagrok-filtering skill). Asserts run synchronously right after `final`,
// so without this settle window every filter assertion is a race.
const ASSERT_SETTLE_MS = 400;
// A judge that merely didn't refuse is not a pass; require real confidence in the answer.
const JUDGE_SCORE_FLOOR = 0.7;
// An on-demand container can take a minute or two to pull, start and index help before it answers.
const READY_ATTEMPTS = 10;
const READY_RETRY_MS = 15000;
const BENCH_DIR = 'System:AppData/Grokky/benchmarks';
const delay = (ms: number): Promise<void> => new Promise((r) => window.setTimeout(r, ms));
const sanitizeLabel = (label: string): string => (label || 'run').replace(/[^\w.-]/g, '_');
const benchPath = (label: string): string => `${BENCH_DIR}/benchmark-${sanitizeLabel(label)}.json`;
const openTableViews = (): DG.TableView[] => Array.from(grok.shell.tableViews);

// ---- Suite + demo tables ---------------------------------------------------------------

async function loadSuite(): Promise<SuitePrompt[]> {
  const text = await _package.files.readAsText('benchmark/suite.yaml');
  return (jsyaml.load(text) as SuitePrompt[]).filter((p) => p && p.prompt);
}

/** Narrows the suite to a comma-separated list of categories, difficulties, or prompt substrings.
 * A full arm is 44 prompts x reps x ~20 s, so validating a suite change — or probing one tier —
 * needs a way to run a slice without paying for the whole thing. */
function filterSuite(suite: SuitePrompt[], only?: string): SuitePrompt[] {
  const terms = (only ?? '').split(',').map((s) => s.trim().toLowerCase()).filter(Boolean);
  if (!terms.length) return suite;
  return suite.filter((p) => terms.some((term) =>
    p.category.toLowerCase() === term || (p.difficulty ?? 'standard').toLowerCase() === term ||
    p.prompt.toLowerCase().indexOf(term) >= 0));
}

const csvCache = new Map<string, DG.DataFrame>();

async function freshTable(name: string): Promise<DG.DataFrame> {
  if (name === 'randomWalk')
    return grok.data.demo.randomWalk();
  if (!csvCache.has(name))
    csvCache.set(name, await grok.dapi.files.readCsv(`System:DemoFiles/${name}.csv`));
  return csvCache.get(name)!.clone();
}

// ---- One measured turn — mirrors runClaudeStreaming() but headless -----------------------

async function runTurn(client: ClaudeRuntimeClient, prompt: string, view: DG.ViewBase | null,
  model?: ClaudeModel, systemPromptMode?: string, difficulty?: string): Promise<RepTurn> {
  // A bare turn declares no tools: the readiness probe must not be gated, graded or tool-driven.
  const viewTools = systemPromptMode === 'none' ? {defs: [], runners: new Map()} : viewFunctionTools();
  const ctx = view ? buildWorkspaceContext() : '';
  const full = ctx ? ctx + '\n---\n\n' + prompt : prompt;
  const sid = `bench-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const targetView = view ?? grok.shell.v;

  return new Promise<RepTurn>((resolve) => {
    const subs: {unsubscribe: () => void}[] = [];
    const cleanup = () => subs.forEach((s) => s.unsubscribe());
    const forSession = (subj: any, fn: (e: any) => void) =>
      subs.push(subj.subscribe((e: any) => {if (e.sessionId === sid) fn(e);}));

    let tSend = 0;
    let ttftMs: number | null = null;
    let toolRoundTrips = 0;
    const tools: string[] = [];
    let accumulated = '';
    let done = false;

    const finish = (extra: Partial<RepTurn>) => {
      if (done) return;
      done = true;
      window.clearTimeout(timer);
      cleanup();
      resolve({
        ok: !extra.error, ttftMs, totalMs: tSend ? performance.now() - tSend : null,
        toolRoundTrips, tools, ctxChars: full.length, toolDefCount: viewTools.defs.length,
        server: null, assertPass: null, unverified: false, content: accumulated, ...extra,
      });
    };

    const timer = window.setTimeout(() => {
      client.abort(sid);
      finish({error: 'turn timed out'});
    }, turnTimeout(difficulty));

    forSession(client.onChunk, (e) => {
      if (ttftMs === null) ttftMs = performance.now() - tSend;
      accumulated += e.content;
    });
    forSession(client.onToolActivity, (e) => {
      toolRoundTrips++;
      if (e.name) tools.push(e.name);
    });

    forSession(client.onInputRequest, async (e) => {
      // Faithful round-trips so multi-tool turns don't stall (the runtime blocks on these).
      if (e.toolName === 'datagrok_exec') {
        const {value, error} = await executeSingleBlock(e.input.code ?? '', targetView, 0);
        let result: any = error ?
          {success: false, error: error.error} :
          {success: true, ...(value != null ? {returnValue: value} : {})};
        if (!error && e.input.verify?.assertion) {
          const v = await runVerification(e.input.verify.assertion, targetView);
          result = {...result, verified: {passed: v.passed,
            ...(v.observed !== undefined ? {observed: v.observed} : {}), ...(v.error ? {error: v.error} : {})}};
        }
        client.respondToInput(sid, e.requestId, result);
      } else if (e.toolName === 'datagrok_verify') {
        const {passed, observed, error} = await runVerification(e.input.assertion ?? '', targetView);
        client.respondToInput(sid, e.requestId,
          {passed, ...(observed !== undefined ? {observed} : {}), ...(error ? {error} : {})});
      } else if (e.toolName === 'datagrok_show_entities')
        client.respondToInput(sid, e.requestId, {success: true});
      else {
        const runner = viewTools.runners.get(e.toolName);
        if (runner) {
          try {
            const result = await runner(e.input ?? {});
            client.respondToInput(sid, e.requestId, result === undefined ? {success: true} : result);
          } catch (err: any) {
            client.respondToInput(sid, e.requestId, {success: false, error: err.message});
          }
        } else // AskUserQuestion / unknown — can't answer headlessly; unblock with an empty answer.
          client.respondToInput(sid, e.requestId, {answers: {}});
      }
    });

    forSession(client.onFinal, (e) => {
      const content = e.revision === 'replaced' ? e.content : (accumulated || e.content);
      // A turn that streamed nothing, called nothing and reported no usage did not happen — the
      // runtime was still starting. Scoring that as a model failure is how a cold container turns
      // into a red row that looks like a regression.
      const empty = !content && !tools.length && !e.metrics;
      finish({server: e.metrics ?? null, unverified: !!e.unverified, content,
        ...(empty ? {error: 'empty turn — runtime not ready'} : {})});
    });
    forSession(client.onError, (e) => finish({error: e.message}));

    tSend = performance.now();
    try {
      client.send(sid, full, {
        ...(viewTools.defs.length ? {clientTools: viewTools.defs} : {}),
        ...(model ? {model} : {}),
        ...(systemPromptMode ? {systemPromptMode} : {}),
      });
    } catch (err: any) {
      finish({error: err.message});
    }
  });
}

// ---- Deterministic assert + Haiku judge -------------------------------------------------

/** Runs one assert expression. Returns the verdict plus the error text when the expression
 * itself threw — a broken assert must be distinguishable from a model failure, not silently false.
 * The body is async so assertions may `await`; previously `await` was a syntax error (silent false)
 * and an async IIFE returned a truthy Promise (a guaranteed false positive). */
async function evalAssert(expr: string, view: DG.ViewBase | null, t: DG.DataFrame | null,
  before: any, opened: number, tools: string[], openedViews: DG.TableView[],
): Promise<{pass: boolean, error?: string}> {
  try {
    // eslint-disable-next-line no-new-func
    const fn = new Function('grok', 'DG', 'view', 't', 'before', 'opened', 'tools', 'openedViews',
      `return (async () => (${expr}))();`);
    return {pass: !!await fn(grok, DG, view, t, before, opened, tools, openedViews)};
  } catch (e: any) {
    const error = `assert threw: ${e?.message ?? e}`;
    console.warn('benchmark:', error, expr);
    return {pass: false, error};
  }
}

/** What the workspace looked like when an assert failed. Without this a failing row only says
 * "wrong", and every investigation turns into guessing whether the assert or the model was at
 * fault — which is how a benchmark ends up blaming the model for its own bugs. Kept short: enough
 * to tell a mis-bound viewer from a missing one, not a full dump. */
function describeState(view: DG.ViewBase | null, t: DG.DataFrame | null, before: any): string {
  const parts: string[] = [];
  try {
    const tv = view as DG.TableView | null;
    if (tv?.viewers) {
      const viewers = Array.from(tv.viewers).filter((v) => v.type !== DG.VIEWER.GRID).map((v) => {
        // Only the column-binding props matter here, and only the ones that are set — a full
        // property dump would bury the one field that is usually wrong.
        const bound = ['xColumnName', 'yColumnName', 'valueColumnName', 'splitColumnName',
          'category1ColumnName', 'colorColumnName', 'categoryColumnNames']
          .map((p) => {
            try {
              const val = (v.props as any)[p];
              return val == null || (Array.isArray(val) && !val.length) ? null : `${p}=${val}`;
            } catch {
              return null;
            }
          }).filter(Boolean);
        return `${v.type}(${bound.join(', ') || 'no columns bound'})`;
      });
      parts.push(`viewers: ${viewers.length ? viewers.join(' · ') : 'none beyond the grid'}`);
    }
    if (t) {
      const added = t.columns.names().filter((n: string) => before.cols.indexOf(n) < 0);
      const removed = before.cols.filter((n: string) => !t.columns.contains(n));
      if (added.length) parts.push(`added cols: ${added.join(', ')}`);
      if (removed.length) parts.push(`removed cols: ${removed.join(', ')}`);
      parts.push(`rows ${t.rowCount}/${before.rowCount}, filtered-in ${t.filter.trueCount}, selected ${t.selection.trueCount}`);
      if (tv?.grid) {
        const hidden = t.columns.names().filter((n: string) => {
          try {
            return tv.grid.col(n)?.visible === false;
          } catch {
            return false;
          }
        });
        if (hidden.length) parts.push(`grid-hidden: ${hidden.join(', ')}`);
      }
    }
  } catch (e: any) {
    parts.push(`(state capture failed: ${e?.message ?? e})`);
  }
  return parts.join(' | ');
}

const JUDGE_SCHEMA = {
  type: 'object',
  properties: {pass: {type: 'boolean'}, score: {type: 'number'}, reason: {type: 'string'}},
  required: ['pass', 'score'],
};

async function judge(client: ClaudeRuntimeClient, prompt: string, answer: string, rubric: string): Promise<JudgeResult> {
  const jp = `You are grading an AI assistant's answer against a rubric. Be strict but fair.\n\n` +
    `USER PROMPT:\n${prompt}\n\nASSISTANT ANSWER:\n${answer || '(empty)'}\n\nRUBRIC:\n${rubric}\n\n` +
    `Return pass (boolean), score (0..1), and a one-line reason.`;
  // One retry: a judge call that times out (query() now has a timeout; the runtime watchdog can
  // kill a wedged CLI) must not flip a correct answer to failed on transient grounds.
  for (let attempt = 0; ; attempt++) {
    try {
      const r = await client.query(jp, {model: ClaudeModel.Haiku, outputSchema: JUDGE_SCHEMA, systemPromptMode: 'none'});
      return {pass: !!r.pass, score: typeof r.score === 'number' ? r.score : (r.pass ? 1 : 0), reason: r.reason ?? ''};
    } catch (e: any) {
      if (attempt >= 1)
        return {pass: false, score: 0, reason: `judge failed: ${e.message}`};
    }
  }
}

// ---- Stats ------------------------------------------------------------------------------

function stats(values: (number | null | undefined)[]): Stats {
  const xs = values.filter((v): v is number => typeof v === 'number' && !isNaN(v)).sort((a, b) => a - b);
  const n = xs.length;
  if (!n) return {n: 0, median: null, mean: null, min: null, max: null, p90: null, std: null};
  const at = (q: number) => xs[Math.min(n - 1, Math.max(0, Math.round(q * (n - 1))))];
  const mean = xs.reduce((s, v) => s + v, 0) / n;
  const variance = xs.reduce((s, v) => s + (v - mean) * (v - mean), 0) / n;
  return {n, median: at(0.5), mean, min: xs[0], max: xs[n - 1], p90: at(0.9), std: Math.sqrt(variance)};
}

// ---- Orchestration ----------------------------------------------------------------------

async function runPrompt(client: ClaudeRuntimeClient, p: SuitePrompt, reps: number,
  model?: ClaudeModel): Promise<PromptResult> {
  const turns: RepTurn[] = [];
  for (let i = 0; i < reps; i++) {
    const preViews = openTableViews();
    let view: DG.TableView | null = null;
    try {
      if (p.table) {
        view = grok.shell.addTableView(await freshTable(p.table));
        await delay(120);
      }
      const t = view ? view.dataFrame : null;
      const before = {
        cols: t ? t.columns.names() : [],
        rowCount: t ? t.rowCount : 0,
        viewers: view ? Array.from(view.viewers).length : 0,
        tableViews: openTableViews().length,
      };
      const turn = await runTurn(client, p.prompt, view, model, undefined, p.difficulty);
      if (turn.ok && p.assert) {
        await delay(ASSERT_SETTLE_MS); // let async filter/viewer updates land before reading state
        const after = openTableViews();
        const openedViews = after.filter((v) => preViews.indexOf(v) < 0 && v !== view);
        const {pass, error} = await evalAssert(p.assert, view, t, before,
          after.length - before.tableViews, turn.tools, openedViews);
        turn.assertPass = pass;
        if (error) turn.assertError = error;
        if (!pass) turn.failState = describeState(view, t, before);
      }
      turns.push(turn);
    } catch (e: any) {
      turns.push({ok: false, ttftMs: null, totalMs: null, toolRoundTrips: 0, tools: [], ctxChars: 0,
        toolDefCount: 0, server: null, assertPass: null, unverified: false, content: '', error: e.message});
    } finally {
      // Close the bench table and anything the assistant opened, so reps stay isolated.
      for (const v of openTableViews())
        if (preViews.indexOf(v) < 0) try {v.close();} catch {}
      await delay(60);
    }
  }

  const ok = turns.filter((t) => t.ok);
  const assertTurns = turns.filter((t) => t.assertPass !== null);
  const assertPassRate = assertTurns.length ? assertTurns.filter((t) => t.assertPass).length / assertTurns.length : null;
  // Judge every rep that produced text, not just the last one: latency is measured at n=reps, so
  // scoring at n=1 let the same defect flip pass/fail between runs.
  const judged = p.rubric ?
    (await Promise.all(ok.filter((t) => t.content).map((t) => judge(client, p.prompt, t.content, p.rubric!)))) : [];
  const judgeRes: JudgeResult | null = judged.length ? {
    pass: judged.filter((j) => j.pass && j.score >= JUDGE_SCORE_FLOOR).length > judged.length / 2,
    score: judged.reduce((s, j) => s + j.score, 0) / judged.length,
    reason: judged.map((j) => j.reason).find((r) => r) ?? '',
  } : null;
  // A prompt that failed every rep is a FAILURE, not "unscored". Leaving it null dropped it from
  // the denominator, which is how a 19-prompt run reported "16/18".
  const accuracyPass = ok.length === 0 ? false :
    assertPassRate !== null ? assertPassRate >= 0.5 :
      judgeRes ? judgeRes.pass : null;

  return {
    category: p.category, prompt: p.prompt,
    difficulty: p.difficulty ?? 'standard', model: model ?? 'default',
    reps, okCount: ok.length,
    toolsUsed: Array.from(new Set(turns.flatMap((t) => t.tools))),
    unverifiedRate: ok.length ? ok.filter((t) => t.unverified).length / ok.length : 0,
    ttftMs: stats(ok.map((t) => t.ttftMs)),
    totalMs: stats(ok.map((t) => t.totalMs)),
    toolRoundTrips: stats(ok.map((t) => t.toolRoundTrips)),
    inputTokens: stats(ok.map((t) => t.server?.inputTokens ?? null)),
    outputTokens: stats(ok.map((t) => t.server?.outputTokens ?? null)),
    cacheReadTokens: stats(ok.map((t) => t.server?.cacheReadTokens ?? null)),
    costUsd: stats(ok.map((t) => t.server?.costUsd ?? null)),
    ctxChars: ok[0]?.ctxChars ?? 0,
    toolDefCount: ok[0]?.toolDefCount ?? 0,
    assertPassRate, judge: judgeRes, accuracyPass,
    ...(turns.find((t) => t.failState)?.failState ? {failState: turns.find((t) => t.failState)!.failState} : {}),
    errors: [
      ...turns.filter((t) => t.error).map((t) => t.error!),
      ...turns.filter((t) => t.assertError).map((t) => t.assertError!),
    ],
  };
}

// ---- Report -----------------------------------------------------------------------------

const r1 = (v: number | null) => v == null ? '—' : (Math.round(v * 10) / 10).toString();
const r0 = (v: number | null) => v == null ? '—' : Math.round(v).toString();
const usd = (v: number | null) => v == null ? '—' : '$' + v.toFixed(4);

/** One rollup row over an arbitrary slice of the results — shared by the difficulty and
 * category tables so the two can never drift in how they aggregate. */
function rollupRow(name: string, rs: PromptResult[], hasTokens: boolean): string {
  const sc = rs.filter((r) => r.accuracyPass !== null);
  const passed = sc.filter((r) => r.accuracyPass).length;
  const cells = [name, `${rs.length}`,
    r0(stats(rs.map((r) => r.totalMs.median)).median),
    r1(stats(rs.map((r) => r.toolRoundTrips.median)).median),
    `${passed}/${sc.length}` + (sc.length ? ` (${Math.round(100 * passed / sc.length)}%)` : '')];
  if (hasTokens) cells.push(usd(rs.reduce((s, r) => s + (r.costUsd.median ?? 0), 0)));
  return '| ' + cells.join(' | ') + ' |';
}

function rollupTable(title: string, groups: [string, PromptResult[]][], hasTokens: boolean): string[] {
  const cols = ['', 'n', 'Median ms', 'Tools', 'Accuracy', ...(hasTokens ? ['Cost'] : [])];
  const lines = ['## ' + title, '', '| ' + cols.join(' | ') + ' |', '|' + cols.map(() => '---').join('|') + '|'];
  for (const [name, rs] of groups) lines.push(rollupRow(name, rs, hasTokens));
  lines.push('');
  return lines;
}

function buildMarkdown(label: string, when: string, reps: number, results: PromptResult[],
  hasTokens: boolean, model: string, only?: string): string {
  const scored = results.filter((r) => r.accuracyPass !== null);
  const passed = scored.filter((r) => r.accuracyPass).length;
  const totalMed = stats(results.map((r) => r.totalMs.median));
  const ttftMed = stats(results.map((r) => r.ttftMs.median));
  const costSum = results.reduce((s, r) => s + (r.costUsd.median ?? 0), 0);

  const lines: string[] = [];
  lines.push(`# Grokky benchmark — \`${label}\``);
  lines.push('');
  lines.push(`| | |`);
  lines.push(`|---|---|`);
  lines.push(`| **When** | ${when} |`);
  lines.push(`| **Model** | \`${model}\` |`);
  lines.push(`| **Prompts × reps** | ${results.length} × ${reps} = ${results.length * reps} turns |`);
  if (only) lines.push(`| **Suite filter** | \`${only}\` — ⚠️ partial run, not comparable to a full arm |`);
  lines.push(`| **Accuracy** | **${passed}/${scored.length}**` +
    (scored.length ? ` (${Math.round(100 * passed / scored.length)}%)` : '') + ` |`);
  lines.push(`| **Median turn** | ${r0(totalMed.median)} ms |`);
  lines.push(`| **Median TTFT** | ${r0(ttftMed.median)} ms |`);
  if (hasTokens) {
    lines.push(`| **Cost** | ${usd(costSum)} total` +
      (passed ? ` · ${usd(costSum / passed)} per passed prompt` : '') + ` |`);
  } else
    lines.push(`| **Cost** | ⚠️ unavailable — rebuild the claude-runtime image to forward SDK usage |`);
  lines.push('');

  // Per-difficulty rollup: a blended number hides where a configuration actually differs. This
  // is what showed a cheap tier taking 8 tool round-trips on prompts the strong tier does in 1.
  const tiers = (['trivial', 'standard', 'hard'] as const).filter((d) => results.some((r) => r.difficulty === d));
  if (tiers.length > 1)
    lines.push(...rollupTable('By difficulty', tiers.map((d) => [d, results.filter((r) => r.difficulty === d)]), hasTokens));

  const cats = Array.from(new Set(results.map((r) => r.category))).sort();
  if (cats.length > 1)
    lines.push(...rollupTable('By category', cats.map((c) => [c, results.filter((r) => r.category === c)]), hasTokens));

  lines.push('## Per prompt');
  lines.push('');
  const cols = hasTokens ?
    ['Category', 'Diff', 'Prompt', 'n', 'TTFT ms', 'Total ms', 'In tok', 'Out tok', 'Cache rd', 'Cost', 'Tools', 'Acc'] :
    ['Category', 'Diff', 'Prompt', 'n', 'TTFT ms', 'Total ms', 'Tools', 'Acc'];
  lines.push('| ' + cols.join(' | ') + ' |');
  lines.push('|' + cols.map(() => '---').join('|') + '|');
  for (const r of results) {
    const prompt = r.prompt.length > 48 ? r.prompt.slice(0, 47) + '…' : r.prompt;
    // Show the raw rate, not a bare tick: 2-of-3 and 3-of-3 are different signals, and the
    // threshold alone hid flaky prompts.
    const rate = r.assertPassRate !== null ? ` ${Math.round(100 * r.assertPassRate)}%` :
      r.judge ? ` ${r.judge.score.toFixed(2)}` : '';
    const acc = r.accuracyPass == null ? '—' : `${r.accuracyPass ? '✅' : '❌'}${rate}`;
    const base = [r.category, r.difficulty, prompt, `${r.okCount}/${r.reps}`,
      r0(r.ttftMs.median), r0(r.totalMs.median)];
    const row = hasTokens ?
      [...base, r0(r.inputTokens.median), r0(r.outputTokens.median), r0(r.cacheReadTokens.median),
        usd(r.costUsd.median), r1(r.toolRoundTrips.median), acc] :
      [...base, r1(r.toolRoundTrips.median), acc];
    lines.push('| ' + row.join(' | ') + ' |');
  }
  lines.push('');
  lines.push('_Values are per-prompt medians across reps. TTFT = time to first streamed token; ' +
    'Total = send→final wall-clock; Tools = tool round-trips/turn. Acc shows the assert pass ' +
    'rate across reps, or the mean judge score._');

  // A red cell in the table above says *that* something failed; this says *what*, which is the
  // difference between a report you can act on and one you have to re-run to understand.
  const failed = results.filter((r) => r.accuracyPass === false);
  if (failed.length) {
    lines.push('');
    lines.push(`## ❌ Failures (${failed.length})`);
    lines.push('');
    for (const r of failed) {
      const why = r.okCount === 0 ? `every rep errored — ${r.errors[0] ?? 'no detail'}` :
        r.assertPassRate !== null ? `assert passed ${Math.round(100 * r.assertPassRate)}% of reps` :
          r.judge ? `judge ${r.judge.score.toFixed(2)} — ${r.judge.reason}` : 'unscored';
      lines.push(`- **${r.category}/${r.difficulty}** \`${r.prompt}\``);
      lines.push(`  - ${why}`);
      if (r.failState) lines.push(`  - left: ${r.failState}`);
      if (r.toolsUsed.length) lines.push(`  - tools: ${r.toolsUsed.join(' → ')}`);
      const errs = Array.from(new Set(r.errors)).slice(0, 2);
      for (const e of errs) lines.push(`  - \`${e.slice(0, 200)}\``);
    }
  }

  // Kept separate from failures on purpose: an assert that throws is *my* bug, and counting it as
  // a model failure is how a harness quietly reports progress it never made.
  const broken = results.filter((r) => r.errors.some((e) => e.startsWith('assert threw')));
  if (broken.length) {
    lines.push('');
    lines.push('## ⚠️ Broken assertions (harness bug, not a model failure)');
    lines.push('');
    for (const r of broken)
      lines.push(`- \`${r.prompt.slice(0, 60)}\` — ${r.errors.find((e) => e.startsWith('assert threw'))}`);
  }
  return lines.join('\n');
}

function downloadText(filename: string, text: string): void {
  const blob = new Blob([text], {type: 'application/octet-stream'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  window.setTimeout(() => {document.body.removeChild(a); URL.revokeObjectURL(url);}, 0);
}

/** Runs the benchmark suite and downloads a `benchmark-<label>.{json,md}` report.
 * `model` pins every turn to one tier ('haiku' | 'sonnet' | 'opus') — that is how the control arms
 * for a model comparison are produced without redeploying the runtime. Omit for the runtime default. */
export async function runBenchmark(label: string, reps: number, model?: string,
  only?: string): Promise<string> {
  label = sanitizeLabel(label);
  reps = reps && reps > 0 ? Math.floor(reps) : 3;
  const pinned = model as ClaudeModel | undefined;
  const client = ClaudeRuntimeClient.getInstance();
  await client.ensureConnected();
  const all = await loadSuite();
  const suite = filterSuite(all, only);
  if (!suite.length)
    return only ? `No suite prompt matches "${only}" (${all.length} in suite).` : 'Benchmark suite is empty.';

  grok.shell.info(`Benchmark "${label}": ${suite.length} prompts × ${reps} reps — warming up…`);
  // The container is on-demand: `ensureConnected` returns as soon as the socket is up, but the
  // runtime behind it can still be booting and will answer instantly with nothing. Block until a
  // turn actually produces output, or the whole arm scores against a container that isn't there.
  for (let i = 0; i < READY_ATTEMPTS; i++) {
    // systemPromptMode 'none' is essential, not an optimisation: under the full prompt the
    // GroundingGate blocks a content-free probe and sends the model off reading help pages until
    // the turn timeout, so the probe measures the gate rather than whether the runtime is up.
    const probe = await runTurn(client, 'Reply with exactly: READY', null, ClaudeModel.Haiku, 'none');
    if (probe.ok && probe.content.trim()) break;
    if (i === READY_ATTEMPTS - 1)
      return `Runtime never became ready (${READY_ATTEMPTS} probes) — is the claude-runtime container up?`;
    grok.shell.info(`Benchmark "${label}": runtime not ready, retrying (${i + 1}/${READY_ATTEMPTS})…`);
    await delay(READY_RETRY_MS);
  }
  // Warm-up turn (discarded) — primes the static-prefix prompt cache with the real system prompt.
  try {
    const w = suite[0];
    const preViews = openTableViews();
    const view = w.table ? grok.shell.addTableView(await freshTable(w.table)) : null;
    await runTurn(client, w.prompt, view, pinned);
    for (const v of openTableViews()) if (preViews.indexOf(v) < 0) try {v.close();} catch {}
  } catch (e) {console.warn('benchmark: warm-up failed', e);}

  const results: PromptResult[] = [];
  for (let i = 0; i < suite.length; i++) {
    grok.shell.info(`Benchmark "${label}": ${i + 1}/${suite.length} — ${suite[i].category}`);
    results.push(await runPrompt(client, suite[i], reps, pinned));
  }

  const hasTokens = results.some((r) => r.inputTokens.n > 0);
  const when = new Date().toISOString();
  const report: BenchReport = {label, when, reps, model: pinned ?? 'default', hasTokens,
    ...(only ? {only} : {}), results};
  const markdown = buildMarkdown(label, when, reps, results, hasTokens, pinned ?? 'default', only);
  downloadText(`benchmark-${label}.json`, JSON.stringify(report, null, 2));
  downloadText(`benchmark-${label}.md`, markdown);
  // Persist server-side too, so compareBenchmarks() can diff runs across configs and a headless
  // runner can retrieve the report without touching the browser's download directory.
  try {
    await grok.dapi.files.writeAsText(benchPath(label), JSON.stringify(report));
    await grok.dapi.files.writeAsText(benchPath(label).replace(/\.json$/, '.md'), markdown);
  } catch (e) {console.warn('benchmark: could not persist report to AppData', e);}

  const scored = results.filter((r) => r.accuracyPass !== null);
  const passed = scored.filter((r) => r.accuracyPass).length;
  const msg = `Benchmark "${label}" complete — ${suite.length} prompts × ${reps} reps, ` +
    `accuracy ${passed}/${scored.length}. Report downloaded.` +
    (hasTokens ? '' : ' (Token metrics unavailable — rebuild claude-runtime image.)');
  grok.shell.info(msg);
  return msg;
}

// ---- Compare two runs -------------------------------------------------------------------

async function loadReport(label: string): Promise<BenchReport | null> {
  try {
    const path = benchPath(label);
    if (!await grok.dapi.files.exists(path)) return null;
    return JSON.parse(await grok.dapi.files.readAsText(path)) as BenchReport;
  } catch {
    return null;
  }
}

function overall(rep: BenchReport) {
  const scored = rep.results.filter((r) => r.accuracyPass !== null);
  return {
    totalMs: stats(rep.results.map((r) => r.totalMs.median)).median,
    ttftMs: stats(rep.results.map((r) => r.ttftMs.median)).median,
    inputTokens: stats(rep.results.map((r) => r.inputTokens.median)).median,
    costSum: rep.results.reduce((s, r) => s + (r.costUsd.median ?? 0), 0),
    accPassed: scored.filter((r) => r.accuracyPass).length,
    accScored: scored.length,
  };
}

const pct = (a: number | null, b: number | null): string =>
  a == null || b == null || a === 0 ? '—' : (b >= a ? '+' : '') + Math.round(100 * (b - a) / a) + '%';

/** Compares two or more saved runs (comma-separated labels) into one Markdown report.
 * N-way rather than pairwise because a model-routing question has at least three arms — cheap,
 * strong, and routed — and reading it as two separate diffs hides where routing actually landed.
 * The first label is the reference every change column is measured against. */
export async function compareBenchmarks(labels: string): Promise<string> {
  const wanted = labels.split(',').map((s) => s.trim()).filter(Boolean);
  if (wanted.length < 2) return 'Give at least two comma-separated run labels to compare.';
  const loaded = await Promise.all(wanted.map(loadReport));
  const missing = wanted.filter((_, i) => !loaded[i]);
  if (missing.length) return `Missing report(s): ${missing.map(sanitizeLabel).join(', ')}`;
  const reps = loaded as BenchReport[];
  const [ref] = reps;
  const os = reps.map(overall);

  const lines: string[] = [];
  lines.push(`# Benchmark comparison — ${reps.map((r) => '`' + r.label + '`').join(' vs ')}`);
  lines.push('');
  lines.push(`Reference arm: \`${ref.label}\`. Change columns are relative to it; ` +
    `for latency and cost, negative is better.`);
  lines.push('');

  const head = ['Metric', ...reps.map((r, i) => r.label + (i ? '' : ' (ref)'))];
  lines.push('| ' + head.join(' | ') + ' |');
  lines.push('|' + head.map(() => '---').join('|') + '|');
  const row = (name: string, get: (o: ReturnType<typeof overall>) => number | null, fmt: (v: number | null) => string) =>
    lines.push(`| ${name} | ` + os.map((o, i) =>
      fmt(get(o)) + (i ? ` (${pct(get(os[0]), get(o))})` : '')).join(' | ') + ' |');
  lines.push(`| Model | ` + reps.map((r) => '`' + r.model + '`').join(' | ') + ' |');
  lines.push(`| Accuracy | ` + os.map((o) => `**${o.accPassed}/${o.accScored}**` +
    (o.accScored ? ` (${Math.round(100 * o.accPassed / o.accScored)}%)` : '')).join(' | ') + ' |');
  row('Median turn (ms)', (o) => o.totalMs, r0);
  row('Median TTFT (ms)', (o) => o.ttftMs, r0);
  row('Median input tokens', (o) => o.inputTokens, r0);
  row('Cost (sum of medians)', (o) => o.costSum, usd);
  lines.push('');

  // Where routing actually landed: a cheap arm is only interesting if it holds accuracy on the
  // easy tiers, so the per-tier accuracy split is the decision-relevant view, not the aggregate.
  const tiers = (['trivial', 'standard', 'hard'] as const)
    .filter((d) => ref.results.some((r) => r.difficulty === d));
  if (tiers.length > 1) {
    lines.push('## Accuracy by difficulty');
    lines.push('');
    lines.push('| Difficulty | ' + reps.map((r) => r.label).join(' | ') + ' |');
    lines.push('|' + ['---', ...reps.map(() => '---')].join('|') + '|');
    for (const d of tiers) {
      lines.push(`| ${d} | ` + reps.map((rep) => {
        const sc = rep.results.filter((r) => r.difficulty === d && r.accuracyPass !== null);
        return `${sc.filter((r) => r.accuracyPass).length}/${sc.length}`;
      }).join(' | ') + ' |');
    }
    lines.push('');
  }

  // Per-prompt, showing accuracy alongside time: a faster arm that quietly stopped passing a
  // prompt is the failure mode a latency-only table is blind to.
  lines.push('## Per prompt — median turn ms (accuracy)');
  lines.push('');
  lines.push('| Category | Diff | Prompt | ' + reps.map((r) => r.label).join(' | ') + ' |');
  lines.push('|' + ['---', '---', '---', ...reps.map(() => '---')].join('|') + '|');
  const indexes = reps.map((rep) => new Map(rep.results.map((r) => [r.prompt, r])));
  for (const rr of ref.results) {
    const prompt = rr.prompt.length > 46 ? rr.prompt.slice(0, 45) + '…' : rr.prompt;
    const cells = indexes.map((ix) => {
      const r = ix.get(rr.prompt);
      if (!r) return '—';
      const mark = r.accuracyPass == null ? '' : r.accuracyPass ? ' ✅' : ' ❌';
      return `${r0(r.totalMs.median)}${mark}`;
    });
    lines.push(`| ${rr.category} | ${rr.difficulty} | ${prompt} | ` + cells.join(' | ') + ' |');
  }

  const md = lines.join('\n');
  const name = `benchmark-compare-${wanted.map(sanitizeLabel).join('-vs-')}`;
  downloadText(`${name}.md`, md);
  try {
    await grok.dapi.files.writeAsText(`${BENCH_DIR}/${name}.md`, md);
  } catch (e) {console.warn('benchmark: could not persist comparison to AppData', e);}
  grok.shell.info(`Compared ${wanted.length} runs. Report downloaded.`);
  return md;
}
