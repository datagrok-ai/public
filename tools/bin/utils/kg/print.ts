/// The terminal rendering of answers and reports, `--output table|csv|json|quiet` over the sections of
/// `answer.ts` and `report.ts`. Nothing here shapes data: the browser and an MCP server take the sections as they are.
import {OutputFormat, printOutput} from '../server-output';
import {Answer, EdgeGroup} from './answer';
import {Report, ReportFormat, markdown} from './report';

/** An op prints wider cells than `grok s`: its lists are ids, and half an id is worse than none. */
const CELL_BUDGET = 80;

/** `via` is the chain a reader follows; `path` is the same chain as data. A table shows one, JSON the other. */
function forOutput(rows: Record<string, unknown>[], output: OutputFormat): Record<string, unknown>[] {
  if (!rows.some((r) => r.path !== undefined)) return rows;
  const drop = output === 'json' ? 'via' : 'path';
  return rows.map((row) => {
    const {[drop]: _, ...rest} = row;
    return rest;
  });
}

/** A cell holds one line: an edge group's lists become their summary, `1` or `0.8–1`, the first evidence, the target ids. */
function flat(row: Record<string, unknown>): Record<string, unknown> {
  const g = row as EdgeGroup;
  const confidence = g.confidence === null ? '' : g.confidence[0] === g.confidence[1] ? String(g.confidence[0]) : `${g.confidence[0]}–${g.confidence[1]}`;
  const more = g.count > g.targets.length ? `, +${g.count - g.targets.length} more` : '';
  return {...g, derived_by: g.derived_by.join(', '), confidence, evidence: g.evidence[0] ?? '', targets: g.targets.map((t) => t.id).join(', ') + more};
}

export function printAnswer(answer: Answer, output: OutputFormat): void {
  if (output === 'json') {
    printOutput({...answer, sections: answer.sections.map((s) => ({...s, rows: forOutput(s.rows, output)}))}, 'json');
    return;
  }
  for (const note of answer.notes ?? []) console.log(note);
  if (answer.target && answer.op !== 'find') console.log(`${answer.op} ${answer.target.id}${answer.target.name ? ` (${answer.target.name})` : ''}`);
  for (const section of answer.sections) {
    if (!section.rows.length) {
      console.log(`\n${section.empty ?? `${section.title} (0)`}`);
      continue;
    }
    const total = section.total ?? section.rows.length;
    const count = section.rows.length < total ? `${section.rows.length} of ${total}; --limit to see more` : `${total}`;
    console.log(`\n${section.title} (${count})`);
    const rows = answer.op === 'explain' && section.title === 'edges' ? section.rows.map(flat) : section.rows;
    if (output !== 'table' || rows[0].group === undefined) {
      printOutput(forOutput(rows, output), output, CELL_BUDGET);
      continue;
    }
    // a table gets the group as a sub-header instead of a repeated column; the rows arrive already in group order
    for (const [name, group] of byGroup(rows)) {
      console.log(`\n  ${name || 'ungrouped'}`);
      printOutput(group.map(({group: _, ...rest}) => rest), output, CELL_BUDGET);
    }
  }
}

function byGroup(rows: Record<string, unknown>[]): [string, Record<string, unknown>[]][] {
  const out: [string, Record<string, unknown>[]][] = [];
  for (const row of rows) {
    const name = String(row.group ?? '');
    if (!out.length || out[out.length - 1][0] !== name) out.push([name, []]);
    out[out.length - 1][1].push(row);
  }
  return out;
}

export function printReport(report: Report, format: ReportFormat): void {
  if (format === 'json') {
    printOutput(report, 'json');
    return;
  }
  if (format === 'md') {
    process.stdout.write(markdown(report));
    return;
  }
  console.log(report.summary);
  for (const note of report.notes) console.log(note);
  for (const section of report.sections) {
    console.log(`\n${section.title} (${section.rows.length})`);
    printOutput(section.rows, 'table');
  }
}
