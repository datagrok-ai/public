/* eslint-disable max-len */
// Shared "Conclusion" results column for the group-comparison tests (two-sample t-test, one-way ANOVA).
// See ttest-results-table-spec.md. The verdict (Significant / Not significant) is a test-level property;
// in a multi-row table (Fisher ANOVA) it sits only on the row that carries the test statistic.

import * as DG from 'datagrok-api/dg';

/** Header of the leftmost results column and of its tooltip. */
export const CONCLUSION_COL_NAME = 'Conclusion';
export const CONCLUSION_HEADER_TOOLTIP = 'Null hypothesis testing';

/** Cell labels for the conclusion column. */
export const CONCLUSION_LABEL = {
  significant: 'Significant',
  notSignificant: 'Not significant',
} as const;

/**
 * Conclusion cell colors (Option 1 — neutral accent), applied via categorical column color-coding.
 * The brand color is a deliberate non-judgemental accent: significance carries no "good/bad" meaning,
 * so red/green are avoided.
 */
const CONCLUSION_COLOR = {
  significant: '#534AB7',
  notSignificant: '#94A1B2',
} as const;

const CONCLUSION_COL_WIDTH = 110;

/** The verdict label for the given significance. */
export function conclusionLabel(significant: boolean): string {
  return significant ? CONCLUSION_LABEL.significant : CONCLUSION_LABEL.notSignificant;
}

/**
 * Build the leftmost `Conclusion` string column for a results grid.
 * The column has `rowCount` rows; the verdict sits only at `valueRow` (the row carrying the test
 * statistic), the rest stay empty — a multi-row table has a single test-level verdict.
 */
export function conclusionColumn(significant: boolean, rowCount = 1, valueRow = 0): DG.Column {
  const values = new Array<string>(rowCount).fill('');
  values[valueRow] = conclusionLabel(significant);
  return DG.Column.fromStrings(CONCLUSION_COL_NAME, values);
}

/**
 * Build the leftmost `Conclusion` column with a per-row verdict (Control comparisons).
 * Each row is its own test, so every row carries its own Significant / Not significant label —
 * unlike `conclusionColumn`, where a single test-level verdict sits on one row.
 */
export function conclusionColumnPerRow(significants: readonly boolean[]): DG.Column {
  return DG.Column.fromStrings(CONCLUSION_COL_NAME, significants.map(conclusionLabel));
}

/** Apply text color-coding (Option 1, neutral accent) + bold to the conclusion column of `grid`. */
export function styleConclusionColumn(grid: DG.Grid): void {
  grid.dataFrame.col(CONCLUSION_COL_NAME)!.meta.colors.setCategorical({
    [CONCLUSION_LABEL.significant]: DG.Color.fromHtml(CONCLUSION_COLOR.significant),
    [CONCLUSION_LABEL.notSignificant]: DG.Color.fromHtml(CONCLUSION_COLOR.notSignificant),
  });

  const col = grid.col(CONCLUSION_COL_NAME)!;
  col.isTextColorCoded = true; // color the text, not the cell background
  col.width = CONCLUSION_COL_WIDTH;

  // Bold the verdict text (color stays driven by the color-coding above); skip the empty cells.
  grid.onCellPrepare((cell) => {
    if (cell.isTableCell && cell.tableColumn?.name === CONCLUSION_COL_NAME && cell.value) {
      const base = cell.style.font ? cell.style.font.replace(/^bold\s+/i, '') : '13px Roboto';
      cell.style.font = `bold ${base}`;
    }
  });
}
