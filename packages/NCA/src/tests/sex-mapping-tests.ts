/**
 * Value→sex (M/F/U) mapping — pure-core unit tests (Slice B, AC-B1).
 *
 * @see src/adapter/sex-mapping.ts
 *
 * Covers the honest-uncertainty heuristic (text tokens resolved, ambiguous codes
 * left unassigned, blank → visible U), the whole-column unique scan, the in-place
 * applier (idempotent, INT→STRING), and the identity-map fast-path predicate.
 *
 * **Footgun gate:** `expect(cond === expected, true, msg)` form throughout.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  defaultSexAssignment, uniqueSexValues, applySexMap, sexMapChangesValues,
  type SexCode,
} from '../adapter/sex-mapping';

category('Sex mapping: value→M/F/U core (Slice B)', () => {
  test('AC-B1: defaultSexAssignment resolves unambiguous tokens, leaves codes unassigned', async () => {
    // Male tokens (case-fold + trim).
    for (const t of ['M', 'm', 'male', 'Male', 'MALE', '  male  '])
      expect(defaultSexAssignment(t), 'M', `"${t}" → M`);
    // Female tokens.
    for (const t of ['F', 'f', 'female', 'Female', 'FEMALE'])
      expect(defaultSexAssignment(t), 'F', `"${t}" → F`);
    // Unknown tokens.
    for (const t of ['U', 'u', 'unk', 'unknown', 'undifferentiated'])
      expect(defaultSexAssignment(t), 'U', `"${t}" → U`);
    // Blank → visible U stratum (never a silent drop).
    expect(defaultSexAssignment(''), 'U', 'blank → U');
    expect(defaultSexAssignment('   '), 'U', 'whitespace → U');
    // Ambiguous numeric codes → null (never guessed — no 1=male convention).
    for (const t of ['0', '1', '2', '9'])
      expect(defaultSexAssignment(t) === null, true, `code "${t}" → unassigned`);
    // Unrecognized token → null.
    expect(defaultSexAssignment('intact') === null, true, 'unknown token → unassigned');
  });

  test('AC-B1: uniqueSexValues scans the whole column; NA → blank bucket', async () => {
    const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex',
      ['M', 'Male', 'M', 'F', 'Female', 'F']);
    expect(JSON.stringify(uniqueSexValues(col)), JSON.stringify(['M', 'Male', 'F', 'Female']),
      'distinct, insertion-ordered');

    const withNa = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', null, 'F']);
    const u = uniqueSexValues(withNa);
    expect(u.includes(''), true, 'NA cell collapses to the blank bucket');
  });

  test('AC-B1: applySexMap rewrites in place, keeps unmapped values, INT→STRING', async () => {
    // Coded INT 1/2 → mapped M/F (column type promoted to STRING).
    const intCol = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Sex', [1, 2, 1, 2]);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['a', 'b', 'c', 'd']),
      intCol,
    ]);
    const map = new Map<string, SexCode>([['1', 'M'], ['2', 'F']]);
    applySexMap(df, 'Sex', map);
    const out = df.col('Sex')!;
    expect(out.type === DG.COLUMN_TYPE.STRING, true, 'coded column promoted to STRING');
    expect(out.getString(0), 'M', 'code 1 → M');
    expect(out.getString(1), 'F', 'code 2 → F');

    // Idempotency: re-applying the raw-keyed map to the now-M/F column is a no-op
    // (M/F are not keys → left unchanged), the re-Compute-safety guarantee.
    applySexMap(df, 'Sex', map);
    expect(df.col('Sex')!.getString(0), 'M', 'idempotent: M stays M on re-apply');
    expect(df.col('Sex')!.getString(1), 'F', 'idempotent: F stays F on re-apply');
  });

  test('AC-B2: an unassigned (unmapped) code stays a distinct stratum — never silently merged', async () => {
    // The science purpose of the `sexUnassigned` Compute-gate: an ambiguous code
    // the scientist has NOT resolved must never be silently folded into M or F.
    // Defensively, even if the gate were bypassed, applySexMap leaves an unmapped
    // value as its own raw stratum (map.get(raw) ?? raw) — it does not coerce it.
    // Here `2` is left out of the map (the "unassigned" state); it must survive as
    // a distinct value, not collapse into the mapped M stratum.
    const col = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Sex', [1, 2, 1, 2]);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['a', 'b', 'c', 'd']),
      col,
    ]);
    applySexMap(df, 'Sex', new Map<string, SexCode>([['1', 'M']])); // 2 deliberately unassigned
    const out = df.col('Sex')!;
    expect(out.getString(0), 'M', 'assigned code 1 → M');
    expect(out.getString(1), '2', 'unassigned code 2 stays its own raw stratum (no silent merge)');
    const strata = new Set<string>();
    for (let i = 0; i < out.length; i++) strata.add(out.getString(i));
    expect(strata.size, 2, 'two distinct strata (M and unresolved 2) — never collapsed to one');
    // defaultSexAssignment leaves codes unassigned, which is what drives the gate.
    expect(defaultSexAssignment('2') === null, true, 'code 2 has no default → UI leaves it unassigned');
  });

  test('AC-B1: sexMapChangesValues distinguishes identity from real maps', async () => {
    const identity = new Map<string, SexCode>([['M', 'M'], ['F', 'F']]);
    expect(sexMapChangesValues(identity), false, 'M→M / F→F is a no-op (Slice-A fast path)');
    const coded = new Map<string, SexCode>([['1', 'M'], ['2', 'F']]);
    expect(sexMapChangesValues(coded), true, 'coded map changes values');
    const blank = new Map<string, SexCode>([['', 'U']]);
    expect(sexMapChangesValues(blank), true, 'blank→U changes values');
  });
});
