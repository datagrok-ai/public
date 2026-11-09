// node eval/grade-run.mjs <runDir> [label]: scores <runDir>/<task>[/<k>]/<entry> as eval/README.md describes.
import {existsSync, mkdirSync, readdirSync, statSync, writeFileSync} from 'node:fs';
import {join, relative, resolve} from 'node:path';
import {EVAL, assertFreshBuild, gradeAll, tasks, today, u2Build} from './lib/grade.mjs';

const [given, label] = process.argv.slice(2);
const runDir = given === undefined ? undefined : resolve(given);
if (runDir === undefined || !existsSync(runDir)) {
  console.error('usage: node eval/grade-run.mjs <runDir> [label]  (a folder of <task>/<entry> or <task>/<k>/<entry>)');
  process.exit(2);
}
assertFreshBuild();

const all = tasks();
const isDir = (p) => existsSync(p) && statSync(p).isDirectory();
const entry = (dir, task) => [task.meta.entry, 'solution.ts'].map((f) => join(dir, f)).find((f) => existsSync(f));
/** A task's output folder: its id, or the name the runs before the neutral ids used (`04-hidden-column`). */
const folders = (task) => [task.id, `${task.id.slice(1)}-${task.meta.name}`].map((d) => join(runDir, d)).filter(isDir);
const byNumber = (a, b) => a.localeCompare(b, undefined, {numeric: true});

const jobs = [];
for (const task of all) {
  for (const dir of folders(task)) {
    const one = entry(dir, task);
    const samples = one !== undefined ? [one] : readdirSync(dir).sort(byNumber)
      .map((k) => isDir(join(dir, k)) ? entry(join(dir, k), task) : undefined).filter((f) => f !== undefined);
    jobs.push(...samples.map((file) => ({task, file})));
  }
}
const results = await gradeAll(jobs);

const mean = (xs) => xs.length === 0 ? 0 : xs.reduce((a, b) => a + b, 0) / xs.length;
const sd = (xs) => Math.sqrt(mean(xs.map((x) => (x - mean(xs)) ** 2)));
const full = (s) => s !== undefined && s.passed === s.total;

const report = all.map((task) => {
  const samples = results.filter((r) => r.task === task.id).map((r) => ({solution: relative(runDir, r.file),
    passed: r.passed, total: r.total, score: r.score,
    failed: Object.fromEntries(r.checks.filter((c) => !c.pass).map((c) => [c.id, c.detail ?? '']))}));
  const scores = samples.map((s) => s.score);
  return {task: task.id, title: task.meta.title, k: samples.length, score: mean(scores),
    pass1: mean(samples.map((s) => full(s) ? 1 : 0)),
    range: samples.length === 0 ? null : [Math.min(...scores), Math.max(...scores)], samples};
});

const k = Math.max(0, ...report.map((t) => t.k));

/** The set's means, and per sample index i the run that one agent per task would have scored (a missing sample
 * counts 0): its spread is how far one run's number can sit from the mean. */
function summary(set) {
  const rounds = Array.from({length: k}, (_, i) => ({score: mean(set.map((t) => t.samples[i]?.score ?? 0)),
    pass1: mean(set.map((t) => full(t.samples[i]) ? 1 : 0))}));
  const spread = (key) => rounds.length === 0 ? null : {min: Math.min(...rounds.map((r) => r[key])),
    max: Math.max(...rounds.map((r) => r[key])), sd: sd(rounds.map((r) => r[key]))};
  return {tasks: set.length, score: mean(set.map((t) => t.score)), pass1: mean(set.map((t) => t.pass1)),
    fullyPassed: set.filter((t) => t.k > 0 && t.pass1 === 1).length, rounds,
    scoreSpread: spread('score'), pass1Spread: spread('pass1')};
}

const summaries = {
  single: summary(report.filter((t) => t.task.startsWith('t'))),
  feature: summary(report.filter((t) => t.task.startsWith('f'))),
  all: summary(report),
};

for (const t of report) {
  const runs = t.samples.map((s) => `${s.passed}/${s.total}`).join(' ') || 'no output';
  console.log(`${t.task.padEnd(5)} ${t.score.toFixed(2)}  pass@1 ${t.pass1.toFixed(2)}  ${runs}`);
  for (const s of t.samples) {
    for (const [id, detail] of Object.entries(s.failed))
      console.log(`    ${s.solution} ${id}: ${detail}`);
  }
}
const range = (x) => x === null ? '' : `, one run ${x.min.toFixed(2)}–${x.max.toFixed(2)} (sd ${x.sd.toFixed(3)})`;
console.log('');
for (const [name, s] of Object.entries(summaries)) {
  console.log(`${name.padEnd(8)} ${String(s.tasks).padStart(2)} tasks  score ${s.score.toFixed(3)}${range(s.scoreSpread)}` +
    `  pass@1 ${s.pass1.toFixed(3)}${range(s.pass1Spread)}  ${s.fullyPassed} passed in full`);
}
console.log(`(${k} sample${k === 1 ? '' : 's'} per task at most; a task with no output counts 0)`);
const scored = JSON.stringify({date: today(), u2: u2Build(), k, score: summaries.all.score,
  fullyPassed: summaries.all.fullyPassed, summary: summaries, tasks: report}, null, 2) + '\n';
writeFileSync(join(runDir, 'score.json'), scored);
if (label !== undefined) {
  mkdirSync(join(EVAL, 'results'), {recursive: true});
  writeFileSync(join(EVAL, 'results', `${today()}-${label}.json`), scored);
}
