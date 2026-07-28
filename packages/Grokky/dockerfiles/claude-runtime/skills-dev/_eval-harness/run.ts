#!/usr/bin/env node
import * as fs from 'node:fs';
import * as path from 'node:path';
import * as url from 'node:url';
import {buildSystemPrompt} from './src/prompt-assembly.ts';
import {invokeClaude} from './src/invoke-claude.ts';
import {evaluate, type PromptSpec, type RubricResult} from './src/rubric.ts';
import {renderReport, type RunMeta} from './src/report.ts';

const __filename = url.fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

function fail(msg: string, code = 1): never {
  process.stderr.write(`[grokky-eval] ${msg}\n`);
  process.exit(code);
}

function timestampSlug(d = new Date()): string {
  const pad = (n: number) => String(n).padStart(2, '0');
  return `${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())}-${pad(d.getHours())}${pad(d.getMinutes())}`;
}

async function main() {
  const skillArg = process.argv[2];
  if (!skillArg)
    fail('Usage: tsx run.ts <skill-dir>\n  e.g. tsx run.ts ../df-and-columns');

  const skillDir = path.resolve(process.cwd(), skillArg);
  if (!fs.existsSync(skillDir) || !fs.statSync(skillDir).isDirectory())
    fail(`Skill directory not found: ${skillDir}`);

  const skillName = path.basename(skillDir);
  const promptsPath = path.join(skillDir, 'eval/prompts.json');
  if (!fs.existsSync(promptsPath))
    fail(`Prompts file not found: ${promptsPath}`);

  const runtimeRoot = path.resolve(__dirname, '..', '..');
  const startedAt = new Date();
  const runDir = path.join(skillDir, 'eval/runs', timestampSlug(startedAt));
  fs.mkdirSync(runDir, {recursive: true});

  // Auth check. The SDK accepts either `ANTHROPIC_API_KEY` or OAuth credentials
  // stored by the Claude Code CLI (macOS keychain, Linux file). If neither path
  // is plausibly set, abort with a stub report — better to surface this loudly
  // than to let the SDK fail with an opaque "process exited with code 1".
  const hasApiKey = !!process.env['ANTHROPIC_API_KEY'];
  if (!hasApiKey) {
    process.stdout.write('[grokky-eval] ANTHROPIC_API_KEY not set — relying on Claude Code OAuth credentials (will fail if unavailable)\n');
  }

  // Build system prompt.
  const assembled = buildSystemPrompt({runtimeRoot, skillDir});
  process.stdout.write(`[grokky-eval] system prompt assembled: ${assembled.systemPrompt.length} chars, ${assembled.inlinedSkillNames.length} skills (${assembled.inlinedSkillNames.join(', ')})\n`);

  // Load prompts.
  const raw = fs.readFileSync(promptsPath, 'utf8');
  const specs = JSON.parse(raw) as PromptSpec[];
  process.stdout.write(`[grokky-eval] running ${specs.length} prompts against ${skillName}\n`);

  const t0 = Date.now();
  const results: RubricResult[] = [];
  const durations = new Map<string, number>();
  const errors = new Map<string, string | undefined>();
  const specsById = new Map<string, PromptSpec>();
  for (const s of specs) specsById.set(s.id, s);

  // Concurrency: keep modest to avoid rate limits on a fresh key.
  const CONCURRENCY = 4;
  let cursor = 0;
  async function worker() {
    while (cursor < specs.length) {
      const idx = cursor++;
      const spec = specs[idx];
      const tStart = Date.now();
      process.stdout.write(`[grokky-eval] [${idx + 1}/${specs.length}] ${spec.id}...\n`);
      let response;
      try {
        response = await invokeClaude({
          prompt: spec.prompt,
          systemPrompt: assembled.systemPrompt,
          model: 'sonnet',
          timeoutMs: 120_000,
        });
      } catch (e: any) {
        response = {text: '', toolUses: [], durationMs: Date.now() - tStart, error: String(e?.message ?? e)};
      }
      durations.set(spec.id, response.durationMs);
      errors.set(spec.id, response.error);
      const result = evaluate(spec, response.text);
      results.push(result);
      const verdict = result.overall ? 'PASS' : 'FAIL';
      process.stdout.write(`[grokky-eval] [${idx + 1}/${specs.length}] ${spec.id}: ${verdict} (${(response.durationMs / 1000).toFixed(1)}s)\n`);
      // Persist raw response next to the report for forensic review.
      const rawPath = path.join(runDir, 'responses');
      fs.mkdirSync(rawPath, {recursive: true});
      fs.writeFileSync(path.join(rawPath, `${spec.id}.txt`), response.text || (response.error ? '[ERROR] ' + response.error : '[empty response]'));
    }
  }
  await Promise.all(Array.from({length: CONCURRENCY}, () => worker()));

  // Sort results by original spec order for stable report.
  results.sort((a, b) => specs.findIndex((s) => s.id === a.id) - specs.findIndex((s) => s.id === b.id));

  const totalDurationMs = Date.now() - t0;
  const meta: RunMeta = {
    skillName,
    skillDir,
    inlinedSkillNames: assembled.inlinedSkillNames,
    model: 'sonnet',
    systemPromptChars: assembled.systemPrompt.length,
    totalDurationMs,
    startedAt: startedAt.toISOString(),
  };

  const report = renderReport({meta, results, specs: specsById, durations, errors});
  const reportPath = path.join(runDir, 'report.md');
  fs.writeFileSync(reportPath, report);

  // Also write a structured JSON next to the report for later diffing.
  fs.writeFileSync(
    path.join(runDir, 'results.json'),
    JSON.stringify({meta, results: results.map((r) => ({
      id: r.id, category: r.category, overall: r.overall,
      path: r.path, helpers: r.helpers, symbols: r.symbols, forbidden: r.forbidden,
      execBlockCount: r.execBlocks.length,
    }))}, null, 2),
  );

  const passed = results.filter((r) => r.overall).length;
  process.stdout.write(`\n[grokky-eval] done: ${passed}/${results.length} passed in ${(totalDurationMs / 1000).toFixed(1)}s\n`);
  process.stdout.write(`[grokky-eval] report: ${reportPath}\n`);
}

main().catch((e) => {
  process.stderr.write(`[grokky-eval] uncaught: ${e?.stack ?? e}\n`);
  process.exit(1);
});
