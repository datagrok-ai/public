// Compares a parsed table with the reference output JMPReader.jl wrote for it (tools/jmp2json.jl).
// No platform imports, so it also runs outside the browser (the Node pass, tools/jmp-corpus.ts).

import {jmpSecondsToDate} from '../formats/jmp/jmp-parser';
import {JmpColumnKind, JmpTable, jmpValue} from '../formats/jmp/jmp-types';

export interface GoldenColumn {
  name: string;
  type: string;
  values: (number | string | null)[];
}

export interface Golden {
  nrows: number;
  columns: GoldenColumn[];
}

/** Which of our column kinds each JMPReader.jl element type may come back as. */
const KINDS: {[type: string]: JmpColumnKind[]} = {
  'Int8': ['int'], 'Int16': ['int'], 'Int32': ['int'], 'Int64': ['bigint'],
  'Float64': ['float', 'currency', 'latitude', 'longitude'],
  'String': ['string'],
  'Date': ['date'],
  'DateTime': ['datetime'],
  'Millisecond': ['duration'],
  'JMPReader.Rowstate': ['rowstate'],
};

/** Returns one line per mismatch; an empty list means the table matches the reference. */
export function compareToGolden(table: JmpTable, golden: Golden, maxPerColumn = 3): string[] {
  const out: string[] = [];
  if (table.info.nrows !== golden.nrows)
    out.push(`rows: ${table.info.nrows}, expected ${golden.nrows}`);
  const names = table.columns.map((c) => c.info.name);
  const expectedNames = golden.columns.map((c) => c.name);
  if (names.join('\u0000') !== expectedNames.join('\u0000'))
    out.push(`columns: [${names.join(', ')}], expected [${expectedNames.join(', ')}]`);

  for (const g of golden.columns) {
    const col = table.columns.find((c) => c.info.name === g.name);
    if (!col)
      continue;
    if (col.decodeError) {
      out.push(`${g.name}: not decoded (${col.decodeError})`);
      continue;
    }
    if (!(KINDS[g.type] ?? []).includes(col.kind)) {
      out.push(`${g.name}: kind ${col.kind}, expected one of ${KINDS[g.type]?.join('/') ?? '?'} for ${g.type}`);
      continue;
    }
    let bad = 0;
    for (let i = 0; i < g.values.length && bad < maxPerColumn; i++) {
      const actual = render(col.kind, jmpValue(col, i));
      if (actual !== g.values[i]) {
        out.push(`${g.name}[${i}]: ${JSON.stringify(actual)}, expected ${JSON.stringify(g.values[i])}`);
        bad++;
      }
    }
  }
  return out;
}

/** Puts a value in the form tools/jmp2json.jl writes for its type. */
function render(kind: JmpColumnKind, v: number | string | null): number | string | null {
  if (v === null)
    return null;
  switch (kind) {
  case 'datetime': return jmpSecondsToDate(v as number).toISOString().replace(/Z$/, '');
  case 'date': return jmpSecondsToDate(v as number).toISOString().slice(0, 10);
  case 'duration': return Math.round((v as number) * 1000) / 1000;
  default: return v;
  }
}
