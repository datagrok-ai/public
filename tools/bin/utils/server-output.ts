/// Docs: [Grok Dapi](/docs/plans/grok-dapi/)
import type {NodeApiError, BatchResponse} from './node-dapi';

export type OutputFormat = 'table' | 'json' | 'csv' | 'quiet';

export function printOutput(data: any, format: OutputFormat): void {
  if (data === null || data === undefined) {
    if (format !== 'quiet') console.log('(empty)');
    return;
  }
  const rows = Array.isArray(data) ? data : [data];
  switch (format) {
    case 'json':
      console.log(JSON.stringify(data, null, 2));
      break;
    case 'csv':
      printCsv(rows);
      break;
    case 'quiet':
      for (const row of rows) {
        if (typeof row === 'object' && row !== null)
          console.log(row.id ?? row.name ?? JSON.stringify(row));
        else
          console.log(row);
      }
      break;
    default:
      printTable(rows);
  }
}

export function cellStr(v: any): string {
  if (v === null || v === undefined) return '';
  if (typeof v === 'object') {
    if (v.name) return v.name;
    if (v.id) return v.id;
    return JSON.stringify(v).slice(0, 40);
  }
  return String(v);
}

function printTable(rows: any[]): void {
  if (!rows.length) { console.log('(no results)'); return; }
  if (typeof rows[0] !== 'object' || rows[0] === null) {
    for (const r of rows) console.log(r);
    return;
  }
  const keys = getKeys(rows);
  const widths: Record<string, number> = {};
  for (const k of keys) widths[k] = k.length;
  for (const row of rows)
    for (const k of keys)
      widths[k] = Math.max(widths[k], cellStr(row[k]).length);

  const header = keys.map((k) => k.padEnd(widths[k])).join('  ');
  const sep = keys.map((k) => '-'.repeat(widths[k])).join('  ');
  console.log(header);
  console.log(sep);
  for (const row of rows)
    console.log(keys.map((k) => cellStr(row[k]).padEnd(widths[k])).join('  '));
}

function printCsv(rows: any[]): void {
  if (!rows.length) return;
  if (typeof rows[0] !== 'object' || rows[0] === null) {
    for (const r of rows) console.log(csvCell(String(r)));
    return;
  }
  const keys = getKeys(rows);
  console.log(keys.map(csvCell).join(','));
  for (const row of rows)
    console.log(keys.map((k) => csvCell(cellStr(row[k]))).join(','));
}

export function csvCell(s: string): string {
  return s.includes(',') || s.includes('"') || s.includes('\n') ? `"${s.replace(/"/g, '""')}"` : s;
}

export function getKeys(rows: any[]): string[] {
  const seen = new Set<string>();
  const keys: string[] = [];
  for (const row of rows)
    for (const k of Object.keys(row))
      if (!seen.has(k)) { seen.add(k); keys.push(k); }
  return keys.slice(0, 12);
}

export function printBatchOutput(response: BatchResponse, format: OutputFormat): void {
  if (format === 'json') {
    console.log(JSON.stringify(response, null, 2));
    return;
  }
  if (format === 'quiet') {
    const s = response.summary;
    console.log(`${s.total} total: ${s.succeeded} succeeded, ${s.failed} failed, ${s.partial} partial, ${s.skipped} skipped`);
    return;
  }

  // table and csv: print summary then per-op rows
  const s = response.summary;
  console.log(`Summary: ${s.total} total  ${s.succeeded} succeeded  ${s.failed} failed  ${s.partial} partial  ${s.skipped} skipped`);
  if (!response.results.length) return;
  console.log('');

  const rows = response.results.map((r) => {
    const brief: string = r.status === 'error'
      ? (r.error?.error ?? 'error')
      : r.status === 'skipped'
        ? (r.reason ?? 'skipped')
        : r.status === 'partial'
          ? `${r.summary?.succeeded}/${r.summary?.total} succeeded`
          : '';
    return {id: r.id ?? '', action: r.action, status: r.status, detail: brief};
  });

  if (format === 'csv')
    printCsv(rows);
  else
    printTable(rows);

  // Write per-op errors to stderr
  const errors = response.results.filter((r) => r.status === 'error');
  if (errors.length)
    process.stderr.write(JSON.stringify(errors.map((r) => ({id: r.id, action: r.action, error: r.error})), null, 2) + '\n');
}

export function printError(err: any): void {
  const apiErr: NodeApiError | undefined = err?.apiError;
  const out: any = apiErr ?? {error: String(err?.message ?? err)};
  process.stderr.write(JSON.stringify(out, null, 2) + '\n');
}
