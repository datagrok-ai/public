import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {ITEM_TYPE} from './types';

/** Regex that matches `${...}` column refs, accounting for escaped braces `\{`/`\}` inside. */
const COL_REF_REGEX = /\$\{((?:\\.|[^}\\])*)\}/g;

/** Matches the leading `${...}` column ref of a formula (the part that names its axis column). */
const LEADING_COL_REF = /^\s*\$\{(?:\\.|[^}\\])*\}/;

/** Escape a column name so it can be safely embedded in `${...}` — wraps nested braces. */
export const escapeColName = (name: string): string => grok.functions.handleOuterBracketsInColName(name, true);

/** Wrap a column name as a `${...}` reference, escaping nested braces. */
export const wrapCol = (name: string): string => `\${${escapeColName(name)}}`;

/** True if `name` is the axis token `x`/`y` — "whatever is on that axis" rather than a column (GROK-17068). */
export const isAxisToken = (name: string | undefined | null): boolean => /^[xy]$/i.test((name ?? '').trim());

/** Replaces the leading `${...}` column ref of `formula` with a ref to `colName`. */
export const replaceLeadingCol = (formula: string, colName: string): string =>
  formula.replace(LEADING_COL_REF, () => wrapCol(colName));

/** The formula with its `${...}` column refs replaced by plain column names; the default title of a line. */
export function getTitleFromFormula(formula: string): string {
  let title = formula;
  if (title) {
    const matches: string[][] = [...title.matchAll(COL_REF_REGEX)]
      .map((m) => [m[0], grok.functions.handleOuterBracketsInColName(m[1], false)]);
    for (const match of matches)
      title = title!.replace(match[0], match[1]);
  }
  return title;
}

/** Returns what is wrong with the structure of a band formula, or an empty string. */
export function validateBandFormula(value: string): string {
  const bandFormulaHelp = 'Band formula should be in format: ${column} in(min, max), or X/Y in(min, max)';
  const v = value ?? '';
  // The left side is a ${...} column ref or a bare axis token; a token carries no ref, hence `<= 1`.
  const pattern = /^\s*(?:\$\{(?:\\.|[^}\\])*\}|[xXyY])\s+in\s*\(\s*[^,()]+\s*,\s*[^,()]+\s*\)\s*$/;
  const refCount = v.match(COL_REF_REGEX)?.length ?? 0;
  return refCount <= 1 && pattern.test(v) ? '' : bandFormulaHelp;
}

/** Returns what is wrong with the structure of a line formula, or an empty string.
 *  [singleAxisColumnName] is the value column of a single-axis viewer (box plot, histogram, bar chart). */
export function validateLineFormula(value: string, singleAxisColumnName?: string): string {
  // Replace each ${...} col-ref with a sentinel to check structure without
  // being confused by `=` chars inside nested column names (e.g. `${${AGE} = 10}`).
  const sentinel = '\u0001COL\u0001';
  const skeleton = value.replace(COL_REF_REGEX, sentinel);
  if (skeleton.split('=').length !== 2)
    return 'Line formula should be in format: ${x or y column} = expression';

  // LHS must be exactly one ${...}, or a bare axis token.
  const lhsSkeleton = skeleton.split('=')[0].trim();
  const lhsIsBareToken = isAxisToken(lhsSkeleton);
  if (lhsSkeleton !== sentinel && !lhsIsBareToken)
    return 'Left side must be a single column in format ${column}, or the axis token X or Y';

  // Extract column refs from the original value; unescape each to get real names.
  const refs = [...value.matchAll(COL_REF_REGEX)]
    .map((m) => grok.functions.handleOuterBracketsInColName(m[1].trim(), false));
  const lhsColumn = lhsIsBareToken ? lhsSkeleton : refs[0];
  const rhsColumns = lhsIsBareToken ? refs : refs.slice(1);

  // Single-axis viewers only support `${valueColumn} = const` formulas — the LHS must reference the value
  // column, and the RHS must be a constant (no column references or operations between them).
  if (singleAxisColumnName) {
    if (!isAxisToken(lhsColumn) && lhsColumn !== singleAxisColumnName)
      return `Left side must reference the value column \${${singleAxisColumnName}}`;
    if (rhsColumns.length > 0)
      return 'Right side must be a constant value — column references are not allowed on single-axis viewers';
  }

  // If RHS references the same column → invalid. Constant-only RHS (zero column refs)
  // is valid — single-axis viewers only produce those.
  if (rhsColumns.length > 0 && rhsColumns.includes(lhsColumn))
    return 'Line formula should contain different columns on both sides of the "=" sign';

  return '';
}

/** The text of an annotation region shown in the "Formula" column of the items table. */
export function formatAreaFormula(item: DG.AnnotationRegion): string {
  if (!item)
    return '';

  if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
    const region = item as DG.AreaAnnotationRegion;
    return `(${region.x}, ${region.y}): ${JSON.stringify(region.area)}`;
  }

  const region = item as DG.FormulaAnnotationRegion;
  return `${region.formula1}; ${region.formula2}`;
}

/** Writes a text in the {@link formatAreaFormula} format back to the region. */
export function parseAreaFormula(item: DG.AnnotationRegion, text: string): void {
  if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
    const region = item as DG.AreaAnnotationRegion;
    const data = text.replace(`(${region.x}, ${region.y}): `, '');
    try {
      region.area = JSON.parse(data || '[]');
    } catch {
      region.area = [];
    }
  } else {
    const region = item as DG.FormulaAnnotationRegion;
    const data = text.split('; ');
    region.formula1 = data[0] ?? '';
    region.formula2 = data[1] ?? '';
  }
}

/** The column on the [axis] of an area region: its `x` / `y` holds the name followed by the axis map. */
export function getAreaColumn(region: DG.AreaAnnotationRegion, axis: 'x' | 'y'): string {
  const map = region[`${axis}Map`];
  const full = region[axis] ?? '';
  return map && full.endsWith(map) ? full.substring(0, full.length - map.length - 1) : full;
}

export function setAreaColumn(region: DG.AreaAnnotationRegion, axis: 'x' | 'y', colName?: string): void {
  const map = region[`${axis}Map`];
  region[axis] = map && colName ? `${colName} ${map}` : colName;
}
