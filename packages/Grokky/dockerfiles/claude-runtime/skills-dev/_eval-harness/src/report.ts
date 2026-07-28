import type {RubricResult, PromptSpec} from './rubric.ts';

export interface RunMeta {
  skillName: string;
  skillDir: string;
  inlinedSkillNames: string[];
  model: string;
  systemPromptChars: number;
  totalDurationMs: number;
  startedAt: string;
}

function dimSummary(d: {passed: boolean; details: {pattern: string; matched: boolean}[]}, isForbidden = false): string {
  if (d.passed) return 'PASS';
  const failed = d.details.filter((x) => (isForbidden ? x.matched : !x.matched));
  return `FAIL (${failed.map((f) => '`' + f.pattern + '`').join(', ')})`;
}

function categorize(results: RubricResult[]): Map<string, {passed: number; total: number}> {
  const m = new Map<string, {passed: number; total: number}>();
  for (const r of results) {
    const e = m.get(r.category) ?? {passed: 0, total: 0};
    e.total += 1;
    if (r.overall) e.passed += 1;
    m.set(r.category, e);
  }
  return m;
}

export function renderReport(opts: {
  meta: RunMeta;
  results: RubricResult[];
  specs: Map<string, PromptSpec>;
  durations: Map<string, number>;
  errors: Map<string, string | undefined>;
}): string {
  const {meta, results, specs, durations, errors} = opts;
  const totalRun = results.length;
  const passed = results.filter((r) => r.overall).length;
  const passRate = totalRun > 0 ? ((passed / totalRun) * 100).toFixed(1) + '%' : 'n/a';

  const byCat = categorize(results);

  const lines: string[] = [];
  lines.push(`# Eval report — ${meta.skillName}`);
  lines.push('');
  lines.push(`Started: ${meta.startedAt}`);
  lines.push(`Duration: ${(meta.totalDurationMs / 1000).toFixed(1)}s`);
  lines.push(`Model: \`${meta.model}\``);
  lines.push(`System prompt size: ${meta.systemPromptChars.toLocaleString()} chars`);
  lines.push(`Inlined skills: ${meta.inlinedSkillNames.map((n) => '`' + n + '`').join(', ')}`);
  lines.push('');
  lines.push('## Summary');
  lines.push('');
  lines.push(`- **Total prompts:** ${totalRun}`);
  lines.push(`- **Passed:** ${passed}`);
  lines.push(`- **Failed:** ${totalRun - passed}`);
  lines.push(`- **Pass rate:** ${passRate}`);
  lines.push('');
  lines.push('### By category');
  lines.push('');
  lines.push('| Category | Passed | Total | Rate |');
  lines.push('|----------|--------|-------|------|');
  for (const [cat, s] of [...byCat.entries()].sort()) {
    const rate = s.total > 0 ? ((s.passed / s.total) * 100).toFixed(0) + '%' : 'n/a';
    lines.push(`| ${cat} | ${s.passed} | ${s.total} | ${rate} |`);
  }
  lines.push('');

  const lats = results.map((r) => durations.get(r.id) ?? 0).filter((n) => n > 0).sort((a, b) => a - b);
  if (lats.length > 0) {
    const mean = lats.reduce((a, b) => a + b, 0) / lats.length;
    const median = lats[Math.floor(lats.length / 2)];
    const p95 = lats[Math.min(lats.length - 1, Math.floor(lats.length * 0.95))];
    const slowest = [...results]
      .map((r) => ({id: r.id, lat: durations.get(r.id) ?? 0}))
      .sort((a, b) => b.lat - a.lat)
      .slice(0, 3);
    lines.push('### Latency');
    lines.push('');
    lines.push(`- Total: ${(meta.totalDurationMs / 1000).toFixed(1)}s wall`);
    lines.push(`- Mean: ${(mean / 1000).toFixed(1)}s · Median: ${(median / 1000).toFixed(1)}s · p95: ${(p95 / 1000).toFixed(1)}s`);
    lines.push(`- Slowest 3: ${slowest.map((s) => `\`${s.id}\` (${(s.lat / 1000).toFixed(1)}s)`).join(', ')}`);
    lines.push('');
  }

  lines.push('## Configuration notes');
  lines.push('');
  lines.push('- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.');
  lines.push('- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`' + meta.skillName + '`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.');
  lines.push('- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.');
  lines.push('- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.');
  lines.push('');

  lines.push('## Per-prompt results');
  lines.push('');
  lines.push('| id | category | path | helpers | symbols | forbidden | overall | latency |');
  lines.push('|----|----------|------|---------|---------|-----------|---------|---------|');
  for (const r of results) {
    const lat = durations.get(r.id);
    lines.push(`| ${r.id} | ${r.category} | ${r.path.passed ? 'pass' : 'fail'} | ${r.helpers.passed ? 'pass' : 'fail'} | ${r.symbols.passed ? 'pass' : 'fail'} | ${r.forbidden.passed ? 'pass' : 'fail'} | **${r.overall ? 'PASS' : 'FAIL'}** | ${lat ? (lat / 1000).toFixed(1) + 's' : '-'} |`);
  }
  lines.push('');

  // Failure deep dives.
  const failures = results.filter((r) => !r.overall);
  if (failures.length > 0) {
    lines.push('## Failed prompts (deep)');
    lines.push('');
    for (const r of failures) {
      const spec = specs.get(r.id);
      lines.push(`### ${r.id} — ${r.category}`);
      lines.push('');
      lines.push(`**Prompt:** ${spec?.prompt ?? '(unknown)'}`);
      if (spec?.intent) lines.push(`**Intent:** ${spec.intent}`);
      if (spec?.trap) lines.push(`**Trap:** ${spec.trap}`);
      lines.push('');
      lines.push('**Rubric:**');
      lines.push('');
      lines.push(`- path: ${dimSummary(r.path)}`);
      lines.push(`- helpers: ${dimSummary(r.helpers)}`);
      lines.push(`- symbols: ${dimSummary(r.symbols)}`);
      lines.push(`- forbidden: ${dimSummary(r.forbidden, true)}`);
      lines.push('');

      // Hypothesis based on which dim failed.
      const hyp: string[] = [];
      if (!r.path.passed) hyp.push('Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class.');
      if (!r.helpers.passed) {
        const miss = r.helpers.details.filter((d) => !d.matched).map((d) => d.pattern);
        hyp.push(`Claude bypassed the expected helper(s) (${miss.map((p) => '`' + p + '`').join(', ')}). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.`);
      }
      if (!r.symbols.passed) {
        const miss = r.symbols.details.filter((d) => !d.matched).map((d) => d.pattern);
        hyp.push(`Expected symbol(s) absent (${miss.map((p) => '`' + p + '`').join(', ')}). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.`);
      }
      if (!r.forbidden.passed) {
        const hits = r.forbidden.details.filter((d) => d.matched).map((d) => d.pattern);
        hyp.push(`Forbidden pattern(s) matched (${hits.map((p) => '`' + p + '`').join(', ')}). Skill's anti-pattern callout isn't strong enough — consider making the wrong path throw, or surface a more obvious "DO NOT" callout near the relevant example.`);
      }
      lines.push('**Hypothesis:** ' + (hyp.join(' ') || '(no rubric dim failed — bug in scoring)'));
      lines.push('');

      const err = errors.get(r.id);
      if (err) {
        lines.push('**SDK error:** ' + err);
        lines.push('');
      }

      if (r.execBlocks.length > 0) {
        lines.push('**Emitted exec blocks:**');
        for (const b of r.execBlocks) {
          lines.push('');
          lines.push('```javascript');
          lines.push(b);
          lines.push('```');
        }
      } else {
        lines.push('**Emitted exec blocks:** (none)');
      }
      lines.push('');
      lines.push('<details><summary>Full response text</summary>');
      lines.push('');
      lines.push('```');
      lines.push(r.responseText);
      lines.push('```');
      lines.push('');
      lines.push('</details>');
      lines.push('');
    }
  }

  // Suggested improvements — aggregate by repeated failure modes.
  lines.push('## Suggested skill improvements');
  lines.push('');
  const tips: string[] = [];
  const helperMisses = new Map<string, string[]>(); // pattern → list of prompt ids
  const forbiddenHits = new Map<string, string[]>();
  for (const r of failures) {
    for (const d of r.helpers.details) {
      if (!d.matched) {
        const arr = helperMisses.get(d.pattern) ?? [];
        arr.push(r.id);
        helperMisses.set(d.pattern, arr);
      }
    }
    for (const d of r.forbidden.details) {
      if (d.matched) {
        const arr = forbiddenHits.get(d.pattern) ?? [];
        arr.push(r.id);
        forbiddenHits.set(d.pattern, arr);
      }
    }
  }
  for (const [pat, ids] of helperMisses)
    tips.push(`Expected helper \`${pat}\` not reached for: ${ids.join(', ')}. Consider tightening the skill's routing table or example density for this case.`);
  for (const [pat, ids] of forbiddenHits)
    tips.push(`Forbidden pattern \`${pat}\` slipped through for: ${ids.join(', ')}. Anti-pattern callout needs to be more prominent — consider a "WRONG / RIGHT" pair right before the relevant example.`);
  if (tips.length === 0) tips.push('No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.');
  for (const t of tips) lines.push(`- ${t}`);
  lines.push('');

  // Anomalies.
  lines.push('## Anomalies');
  lines.push('');
  const anomalies: string[] = [];
  for (const r of results) {
    const err = errors.get(r.id);
    if (err) anomalies.push(`\`${r.id}\`: SDK error — ${err}`);
    if (r.overall && r.execBlocks.length === 0) anomalies.push(`\`${r.id}\`: passed rubric but emitted no exec blocks (expected_path may not have required one).`);
    if (r.execBlocks.length > 3) anomalies.push(`\`${r.id}\`: emitted ${r.execBlocks.length} exec blocks — unusually verbose; check whether the skill should encourage consolidation.`);
  }
  if (anomalies.length === 0) anomalies.push('None.');
  for (const a of anomalies) lines.push(`- ${a}`);
  lines.push('');

  return lines.join('\n');
}
