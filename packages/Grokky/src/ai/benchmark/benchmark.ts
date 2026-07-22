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
}

interface RepTurn {
  ok: boolean;
  ttftMs: number | null;
  totalMs: number | null;
  toolRoundTrips: number;
  ctxChars: number;
  toolDefCount: number;
  server: TurnServerMetrics | null;
  assertPass: boolean | null;
  content: string;
  error?: string;
}

interface Stats {n: number; median: number | null; mean: number | null; min: number | null; max: number | null; p90: number | null; std: number | null}
interface JudgeResult {pass: boolean; score: number; reason: string}

interface PromptResult {
  category: string;
  prompt: string;
  reps: number;
  okCount: number;
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
  errors: string[];
}

interface BenchReport {label: string; when: string; reps: number; hasTokens: boolean; results: PromptResult[]}

const TURN_TIMEOUT_MS = 120000;
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

const csvCache = new Map<string, DG.DataFrame>();

async function freshTable(name: string): Promise<DG.DataFrame> {
  if (name === 'randomWalk')
    return grok.data.demo.randomWalk();
  if (!csvCache.has(name))
    csvCache.set(name, await grok.dapi.files.readCsv(`System:DemoFiles/${name}.csv`));
  return csvCache.get(name)!.clone();
}

// ---- One measured turn — mirrors runClaudeStreaming() but headless -----------------------

async function runTurn(client: ClaudeRuntimeClient, prompt: string, view: DG.ViewBase | null): Promise<RepTurn> {
  const viewTools = viewFunctionTools();
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
    let accumulated = '';
    let done = false;

    const finish = (extra: Partial<RepTurn>) => {
      if (done) return;
      done = true;
      window.clearTimeout(timer);
      cleanup();
      resolve({
        ok: !extra.error, ttftMs, totalMs: tSend ? performance.now() - tSend : null,
        toolRoundTrips, ctxChars: full.length, toolDefCount: viewTools.defs.length,
        server: null, assertPass: null, content: accumulated, ...extra,
      });
    };

    const timer = window.setTimeout(() => {
      client.abort(sid);
      finish({error: 'turn timed out'});
    }, TURN_TIMEOUT_MS);

    forSession(client.onChunk, (e) => {
      if (ttftMs === null) ttftMs = performance.now() - tSend;
      accumulated += e.content;
    });
    forSession(client.onToolActivity, () => toolRoundTrips++);

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

    forSession(client.onFinal, (e) => finish({server: e.metrics ?? null,
      content: e.revision === 'replaced' ? e.content : (accumulated || e.content)}));
    forSession(client.onError, (e) => finish({error: e.message}));

    tSend = performance.now();
    try {
      client.send(sid, full, {...(viewTools.defs.length ? {clientTools: viewTools.defs} : {})});
    } catch (err: any) {
      finish({error: err.message});
    }
  });
}

// ---- Deterministic assert + Haiku judge -------------------------------------------------

function evalAssert(expr: string, view: DG.ViewBase | null, t: DG.DataFrame | null,
  before: any, opened: number): boolean {
  try {
    // eslint-disable-next-line no-new-func
    const fn = new Function('grok', 'DG', 'view', 't', 'before', 'opened', `return (${expr});`);
    return !!fn(grok, DG, view, t, before, opened);
  } catch (e) {
    console.warn('benchmark: assert threw', expr, e);
    return false;
  }
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
  try {
    const r = await client.query(jp, {model: ClaudeModel.Haiku, outputSchema: JUDGE_SCHEMA, systemPromptMode: 'none'});
    return {pass: !!r.pass, score: typeof r.score === 'number' ? r.score : (r.pass ? 1 : 0), reason: r.reason ?? ''};
  } catch (e: any) {
    return {pass: false, score: 0, reason: `judge failed: ${e.message}`};
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

async function runPrompt(client: ClaudeRuntimeClient, p: SuitePrompt, reps: number): Promise<PromptResult> {
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
        viewers: view ? Array.from(view.viewers).length : 0,
        tableViews: openTableViews().length,
      };
      const turn = await runTurn(client, p.prompt, view);
      if (turn.ok && p.assert) {
        const opened = openTableViews().length - before.tableViews;
        turn.assertPass = evalAssert(p.assert, view, t, before, opened);
      }
      turns.push(turn);
    } catch (e: any) {
      turns.push({ok: false, ttftMs: null, totalMs: null, toolRoundTrips: 0, ctxChars: 0,
        toolDefCount: 0, server: null, assertPass: null, content: '', error: e.message});
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
  const lastContent = [...ok].reverse().find((t) => t.content)?.content ?? '';
  const judgeRes = p.rubric && lastContent ? await judge(client, p.prompt, lastContent, p.rubric) : null;
  const accuracyPass = assertPassRate !== null ? assertPassRate >= 0.5 : (judgeRes ? judgeRes.pass : null);

  return {
    category: p.category, prompt: p.prompt, reps, okCount: ok.length,
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
    errors: turns.filter((t) => t.error).map((t) => t.error!),
  };
}

// ---- Report -----------------------------------------------------------------------------

const r1 = (v: number | null) => v == null ? '—' : (Math.round(v * 10) / 10).toString();
const r0 = (v: number | null) => v == null ? '—' : Math.round(v).toString();
const usd = (v: number | null) => v == null ? '—' : '$' + v.toFixed(4);

function buildMarkdown(label: string, when: string, reps: number, results: PromptResult[], hasTokens: boolean): string {
  const scored = results.filter((r) => r.accuracyPass !== null);
  const passed = scored.filter((r) => r.accuracyPass).length;
  const totalMed = stats(results.map((r) => r.totalMs.median));
  const ttftMed = stats(results.map((r) => r.ttftMs.median));
  const costSum = results.reduce((s, r) => s + (r.costUsd.median ?? 0), 0);

  const lines: string[] = [];
  lines.push(`# Grokky benchmark — \`${label}\``);
  lines.push('');
  lines.push(`- **When:** ${when}`);
  lines.push(`- **Prompts:** ${results.length} · **Reps each:** ${reps}`);
  lines.push(`- **Median turn:** ${r0(totalMed.median)} ms · **Median TTFT:** ${r0(ttftMed.median)} ms`);
  lines.push(`- **Accuracy:** ${passed}/${scored.length} prompts passed` +
    (scored.length ? ` (${Math.round(100 * passed / scored.length)}%)` : ''));
  if (hasTokens) lines.push(`- **Cost (sum of per-prompt medians):** ${usd(costSum)}`);
  else lines.push(`- ⚠️ **Token/cost metrics unavailable** — rebuild the claude-runtime image to forward SDK usage (see docs/BENCHMARK.md).`);
  lines.push('');

  const cols = hasTokens ?
    ['Category', 'Prompt', 'n', 'TTFT ms', 'Total ms', 'In tok', 'Out tok', 'Cache rd', 'Cost', 'Tools', 'Acc'] :
    ['Category', 'Prompt', 'n', 'TTFT ms', 'Total ms', 'Tools', 'Acc'];
  lines.push('| ' + cols.join(' | ') + ' |');
  lines.push('|' + cols.map(() => '---').join('|') + '|');
  for (const r of results) {
    const prompt = r.prompt.length > 48 ? r.prompt.slice(0, 47) + '…' : r.prompt;
    const acc = r.accuracyPass == null ? '—' : (r.accuracyPass ? '✅' : '❌');
    const base = [r.category, prompt, `${r.okCount}/${r.reps}`, r0(r.ttftMs.median), r0(r.totalMs.median)];
    const row = hasTokens ?
      [...base, r0(r.inputTokens.median), r0(r.outputTokens.median), r0(r.cacheReadTokens.median),
        usd(r.costUsd.median), r1(r.toolRoundTrips.median), acc] :
      [...base, r1(r.toolRoundTrips.median), acc];
    lines.push('| ' + row.join(' | ') + ' |');
  }
  lines.push('');
  lines.push('_Values are per-prompt medians across reps. TTFT = time to first streamed token; ' +
    'Total = send→final wall-clock; Tools = tool round-trips/turn._');
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

/** Runs the benchmark suite and downloads a `benchmark-<label>.{json,md}` report. */
export async function runBenchmark(label: string, reps: number): Promise<string> {
  label = sanitizeLabel(label);
  reps = reps && reps > 0 ? Math.floor(reps) : 3;
  const client = ClaudeRuntimeClient.getInstance();
  await client.ensureConnected();
  const suite = await loadSuite();
  if (!suite.length) return 'Benchmark suite is empty.';

  grok.shell.info(`Benchmark "${label}": ${suite.length} prompts × ${reps} reps — warming up…`);
  // Warm-up turn (discarded) — primes the static-prefix prompt cache and the container.
  try {
    const w = suite[0];
    const preViews = openTableViews();
    const view = w.table ? grok.shell.addTableView(await freshTable(w.table)) : null;
    await runTurn(client, w.prompt, view);
    for (const v of openTableViews()) if (preViews.indexOf(v) < 0) try {v.close();} catch {}
  } catch (e) {console.warn('benchmark: warm-up failed', e);}

  const results: PromptResult[] = [];
  for (let i = 0; i < suite.length; i++) {
    grok.shell.info(`Benchmark "${label}": ${i + 1}/${suite.length} — ${suite[i].category}`);
    results.push(await runPrompt(client, suite[i], reps));
  }

  const hasTokens = results.some((r) => r.inputTokens.n > 0);
  const when = new Date().toISOString();
  const report: BenchReport = {label, when, reps, hasTokens, results};
  const markdown = buildMarkdown(label, when, reps, results, hasTokens);
  downloadText(`benchmark-${label}.json`, JSON.stringify(report, null, 2));
  downloadText(`benchmark-${label}.md`, markdown);
  // Persist server-side too, so compareBenchmarks() can diff runs across configs.
  try {
    await grok.dapi.files.writeAsText(benchPath(label), JSON.stringify(report));
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

/** Diffs two saved runs (by label) into a Markdown delta report and downloads it. */
export async function compareBenchmarks(labelA: string, labelB: string): Promise<string> {
  const [a, b] = await Promise.all([loadReport(labelA), loadReport(labelB)]);
  if (!a || !b)
    return `Missing report(s): ${!a ? sanitizeLabel(labelA) : ''} ${!b ? sanitizeLabel(labelB) : ''}`.trim();

  const oa = overall(a);
  const ob = overall(b);
  const lines: string[] = [];
  lines.push(`# Benchmark comparison — \`${a.label}\` → \`${b.label}\``);
  lines.push('');
  lines.push('| Metric | ' + a.label + ' | ' + b.label + ' | Change |');
  lines.push('|---|---|---|---|');
  lines.push(`| Median turn (ms) | ${r0(oa.totalMs)} | ${r0(ob.totalMs)} | ${pct(oa.totalMs, ob.totalMs)} |`);
  lines.push(`| Median TTFT (ms) | ${r0(oa.ttftMs)} | ${r0(ob.ttftMs)} | ${pct(oa.ttftMs, ob.ttftMs)} |`);
  lines.push(`| Median input tokens | ${r0(oa.inputTokens)} | ${r0(ob.inputTokens)} | ${pct(oa.inputTokens, ob.inputTokens)} |`);
  lines.push(`| Cost (sum of medians) | ${usd(oa.costSum)} | ${usd(ob.costSum)} | ${pct(oa.costSum, ob.costSum)} |`);
  lines.push(`| Accuracy | ${oa.accPassed}/${oa.accScored} | ${ob.accPassed}/${ob.accScored} | — |`);
  lines.push('');
  lines.push('## Per-prompt turn time (ms, median)');
  lines.push('');
  lines.push('| Category | Prompt | ' + a.label + ' | ' + b.label + ' | Change |');
  lines.push('|---|---|---|---|---|');
  const byPrompt = new Map(b.results.map((r) => [r.prompt, r]));
  for (const ra of a.results) {
    const rb = byPrompt.get(ra.prompt);
    if (!rb) continue;
    const prompt = ra.prompt.length > 48 ? ra.prompt.slice(0, 47) + '…' : ra.prompt;
    lines.push(`| ${ra.category} | ${prompt} | ${r0(ra.totalMs.median)} | ${r0(rb.totalMs.median)} | ${pct(ra.totalMs.median, rb.totalMs.median)} |`);
  }
  const md = lines.join('\n');
  downloadText(`benchmark-compare-${sanitizeLabel(labelA)}-vs-${sanitizeLabel(labelB)}.md`, md);
  grok.shell.info(`Compared "${a.label}" → "${b.label}". Diff downloaded.`);
  return md;
}
