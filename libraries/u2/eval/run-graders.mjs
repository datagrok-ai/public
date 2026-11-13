// node eval/run-graders.mjs [--repeat N] [task prefix …]: the page keeps what AGENT-PROMPT.md promises; good passes
// all, each bad fails exactly what task.json lists (an `absorbed` one, a bug u2 now contains, fails nothing), in each
// of N rounds (default 1); a verdict that changes between rounds is a flaky check.
import {spawnSync} from 'node:child_process';
import {mkdirSync, readdirSync, writeFileSync} from 'node:fs';
import {join} from 'node:path';
import {EVAL, assertFreshBuild, failedIds, gradeAll, tasks, today, u2Build} from './lib/grade.mjs';

assertFreshBuild();
const probe = spawnSync(process.execPath, ['--import', join(EVAL, 'lib', 'register.mjs'),
  join(EVAL, 'lib', 'promises.mjs')], {encoding: 'utf8', timeout: 60_000});
const unkept = (() => {
  try {
    return JSON.parse(probe.stdout.trim().split('\n').at(-1));
  }
  catch (_) {
    return [`the probe did not finish: ${(probe.stderr || probe.stdout).trim().split('\n').slice(-2).join(' | ')}`];
  }
})();
console.log(unkept.length === 0 ? 'the page keeps every API AGENT-PROMPT.md promises\n' :
  `the page lacks what AGENT-PROMPT.md promises: ${unkept.join('; ')}\n`);
const args = process.argv.slice(2);
const at = args.indexOf('--repeat');
const repeat = at < 0 ? 1 : Number(args.splice(at, 2)[1]);
const only = args;
const chosen = tasks().filter((t) => only.length === 0 || only.some((p) => t.id.startsWith(p)));
const jobs = chosen.flatMap((task) => ['good.ts', ...readdirSync(task.dir).filter((f) => /^bad-.*\.ts$/.test(f)).sort()]
  .map((f) => ({task, file: join(task.dir, f), name: f.slice(0, -3)})));

const reasons = (r) => Object.fromEntries(r.checks.filter((c) => !c.pass).map((c) => [c.id, c.detail ?? '']));

async function round() {
  const results = await gradeAll(jobs);
  const report = chosen.map((task) => {
    const mine = results.map((r, i) => ({...r, name: jobs[i].name})).filter((r) => r.task === task.id);
    const good = mine.find((r) => r.name === 'good');
    const bad = mine.filter((r) => r.name !== 'good').map((r) => {
      const absorbed = task.meta.absorbed?.[r.name] !== undefined;
      const expected = [...task.meta.bad?.[r.name] ?? []].sort();
      const failed = failedIds(r);
      const asExpected = absorbed ? failed.length === 0 : expected.length > 0 && failed.join() === expected.join();
      return {name: r.name, absorbed, expected, failed, asExpected, reasons: reasons(r)};
    });
    const unlisted = [...Object.keys(task.meta.bad ?? {}), ...Object.keys(task.meta.absorbed ?? {})]
      .filter((n) => !bad.some((b) => b.name === n));
    return {task: task.id, title: task.meta.title, source: task.meta.source, checks: good.checks.map((c) => c.id),
      good: {passed: good.passed, total: good.total, failed: failedIds(good), reasons: reasons(good)},
      bad, missingBad: unlisted};
  });
  const badAll = report.flatMap((t) => t.bad);
  const summary = {
    tasks: report.length,
    checks: report.reduce((n, t) => n + t.checks.length, 0),
    goodPassing: report.filter((t) => t.good.failed.length === 0).length,
    badSolutions: badAll.length,
    badFailingAsExpected: badAll.filter((b) => b.asExpected).length,
  };
  const valid = summary.goodPassing === summary.tasks && summary.badFailingAsExpected === summary.badSolutions &&
    report.every((t) => t.missingBad.length === 0);
  return {report, summary, valid, verdicts: results.map(failedIds)};
}

function print({report, summary, valid}, n) {
  for (const t of report) {
    const goodFails = t.good.failed.length ? `  FAILS ${t.good.failed.join(', ')}` : '';
    console.log(`${t.task}  good ${t.good.passed}/${t.good.total}${goodFails}`);
    for (const b of t.bad) {
      const wanted = b.asExpected ? '' : b.absorbed ? ' (listed as absorbed by u2)' : ` (expected ${b.expected.join(', ')})`;
      const said = b.absorbed && b.asExpected ? 'absorbed by u2, fails nothing' : `fails ${b.failed.join(', ') || 'nothing'}`;
      console.log(`  ${b.asExpected ? 'ok  ' : 'BAD '} ${b.name}: ${said}${wanted}`);
    }
    for (const name of t.missingBad)
      console.log(`  BAD  ${name}: listed in task.json, no ${name}.ts`);
  }
  console.log(`\nround ${n}: ${summary.goodPassing}/${summary.tasks} good solutions pass, ` +
    `${summary.badFailingAsExpected}/${summary.badSolutions} bad solutions fail as expected: ` +
    `graders ${valid ? 'valid' : 'NOT valid'}\n`);
}

const rounds = [];
for (let n = 1; n <= repeat; n++) {
  const r = await round();
  rounds.push(r);
  if (n === 1 || !r.valid)
    print(r, n);
  else
    console.log(`round ${n}: graders ${r.valid ? 'valid' : 'NOT valid'}`);
}
const flaky = jobs.map((j, i) => ({job: `${j.task.id}/${j.name}`, seen: [...new Set(rounds.map((r) => r.verdicts[i].join()))]}))
  .filter((f) => f.seen.length > 1);
for (const f of flaky)
  console.log(`FLAKY ${f.job}: failed ${f.seen.map((s) => `[${s}]`).join(' / ')} in different rounds`);
const valid = unkept.length === 0 && rounds.every((r) => r.valid) && flaky.length === 0;
if (repeat > 1)
  console.log(`\n${repeat} rounds: graders ${valid ? 'valid in every round' : 'NOT valid'}`);

const last = rounds.at(-1);
if (only.length === 0) {
  mkdirSync(join(EVAL, 'results'), {recursive: true});
  writeFileSync(join(EVAL, 'results', 'baseline-graders.json'), JSON.stringify({date: today(), u2: u2Build(), valid,
    unkept, rounds: rounds.map((r) => ({valid: r.valid, summary: r.summary})), flaky, summary: last.summary,
    tasks: last.report}, null, 2) + '\n');
}
process.exit(valid ? 0 : 1);
