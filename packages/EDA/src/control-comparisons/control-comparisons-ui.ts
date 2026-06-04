/* eslint-disable max-len */
// Control comparisons (many-to-one) - UI.
//
// Mirrors the t-test / ANOVA dialogs and the shared results pattern: one box plot docked right,
// one results grid docked below (no tab control), a leftmost Conclusion column with per-row
// verdicts (neutral purple/gray accent). Computational core lives in control-comparisons-tools.

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NumCol} from '../anova/anova-tools';
import {
  controlComparisons, ControlComparisonsReport, ControlComparisonsMethod,
} from './control-comparisons-tools';
import {
  CONCLUSION_COL_NAME, CONCLUSION_HEADER_TOOLTIP, conclusionColumnPerRow, styleConclusionColumn,
} from '../group-comparison/conclusion-column';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const HELP_URL = '/help/explore/group-comparison#control-comparisons';

enum SIGNIFICANCE {
  DEFAULT = 0.05,
  MIN = 0.0001,
  MAX = 0.99,
  INFIMUM = 0,
  SUPREMUM = 1,
}

const METHODS: ControlComparisonsMethod[] = ['Dunnett', 'Holm-Welch'];

const METHOD_TOOLTIP = ui.markdown(
  'Set the method for the test:\n\n' +
  '* **Dunnett** — classical procedure for comparing several groups against a single control. ' +
  'Assumes equal variances across groups. More powerful than Holm-Welch because it accounts for ' +
  'the shared control. Standard in toxicology and dose-response studies.\n\n' +
  '* **Holm-Welch** — pairwise Welch\'s t-tests of each group against the control, with Holm ' +
  'correction for multiple comparisons. Robust to unequal variances. Easier to interpret, but ' +
  'slightly less powerful when variances are equal.\n\n',
);

// ── Category helpers ─────────────────────────────────────────────

/** Factorization bin count for a categorical column (upper bound on category codes). */
function binCount(factor: DG.Column): number {
  return factor.type === DG.COLUMN_TYPE.BOOL ? 2 : factor.categories.length;
}

/** Number of non-empty groups (distinct non-null values). */
function groupCount(factor: DG.Column): number {
  return factor.stats.uniqueCount;
}

/** True when every non-null value is distinct (an ID-like column — each group has size 1, so the
 *  column cannot define groups). Mirrors the ANOVA dialog's guard. */
function isFactorAllUnique(factor: DG.Column): boolean {
  const uniqueCount = factor.stats.uniqueCount;
  const nonNullCount = factor.length - factor.stats.missingValueCount;
  return uniqueCount >= 2 && uniqueCount === nonNullCount;
}

const ALL_UNIQUE_MSG =
  'Every value is unique (an ID-like column) — no groups to compare. ' +
  'Pick a column with at least one repeated category.';

/** Human-readable label for a category code. */
function labelOf(factor: DG.Column, code: number): string {
  if (factor.type === DG.COLUMN_TYPE.BOOL)
    return code === 1 ? 'true' : 'false';
  return String(factor.categories[code]);
}

/** Category code for a label. */
function codeOf(factor: DG.Column, label: string): number {
  if (factor.type === DG.COLUMN_TYPE.BOOL)
    return label === 'true' ? 1 : 0;
  return factor.categories.indexOf(label);
}

/** Labels of the (present) groups, in category-code order. */
function groupLabels(factor: DG.Column): string[] {
  if (factor.type === DG.COLUMN_TYPE.BOOL)
    return ['false', 'true'];
  return Array.from(factor.categories);
}

// ── Number formatting ────────────────────────────────────────────

/** Adaptive p-value column format: scientific below 0.001, fixed "0.000" otherwise. */
function pColumnFormat(pValues: readonly number[]): string {
  return pValues.some((p) => p < 0.001) ? 'scientific' : '0.000';
}

// ── Box plot description ─────────────────────────────────────────

/** Box plot caption: "None / ${m} of ${k-1} groups differ significantly from "${control}"". */
function buildDescription(report: ControlComparisonsReport, controlLabel: string): string {
  const total = report.comparisons.length;
  const sig = report.comparisons.filter((c) => c.significant).length;
  // No significant difference: avoid "differ significantly" (not significant ≠ proven equal) —
  // state the absence of a detectable difference instead.
  if (sig === 0)
    return `No detectable difference from "${controlLabel}"`;
  const groups = `group${total === 1 ? '' : 's'}`;
  return `${sig} of ${total} ${groups} differ significantly from "${controlLabel}"`;
}

// ── Per-row Conclusion cell tooltip (H0 / H1 / Conclusion + verdict) ──

/** Null-hypothesis-testing tooltip for a conclusion cell, in the ANOVA dialog's compact style. */
function conclusionTooltip(significant: boolean, groupLabel: string, controlLabel: string,
  featureName: string): HTMLElement {
  return ui.divV([
    ui.markdown(`**H0:** mean "${featureName}" in "${groupLabel}" equals mean in "${controlLabel}".`),
    ui.markdown(`**H1:** mean "${featureName}" in "${groupLabel}" differs from mean in "${controlLabel}".`),
    ui.markdown(`**Conclusion:** ${significant ?
      'Reject the null hypothesis.' :
      'Fail to reject the null hypothesis.'}`),
    ui.markdown(significant ?
      `**"${groupLabel}" differs significantly from "${controlLabel}" in mean "${featureName}".**` :
      `**"${groupLabel}" does not differ significantly from "${controlLabel}" in mean "${featureName}".**`),
  ]);
}

// ── Column header tooltips ───────────────────────────────────────

function headerTooltips(report: ControlComparisonsReport, controlLabel: string, featureName: string,
  ciCaption: string): Map<string, string> {
  const dunnett = report.method === 'Dunnett';
  const diff = `mean(group) − mean("${controlLabel}"), in the original units of "${featureName}"`;

  return new Map<string, string>([
    [CONCLUSION_COL_NAME, CONCLUSION_HEADER_TOOLTIP],
    ['Group', `Group compared against the control "${controlLabel}".`],
    ['n', 'Sample size (n). Number of observations in this group.'],
    ['Mean diff', `Mean difference. ${diff}.`],
    [`${ciCaption} low`, ciTooltip(report.alpha)],
    [`${ciCaption} high`, ciTooltip(report.alpha)],
    ['t', dunnett ?
      'Dunnett\'s t-statistic. Uses the pooled within-group variance; compared against Dunnett\'s ' +
      'multivariate-t reference distribution.' :
      'Welch\'s t-statistic. Uses each group\'s own variance (no pooling).'],
    ['df', dunnett ?
      'Pooled degrees of freedom (N − k). A single value shared across all comparisons.' :
      'Welch–Satterthwaite degrees of freedom. Fractional by design; each comparison has its own.'],
    ['p (raw)', dunnett ?
      'Uncorrected per-comparison p-value. Not yet adjusted for the k−1 comparisons.' :
      'Uncorrected two-sided Welch p-value for this comparison. Not yet adjusted for the k−1 comparisons.'],
    ['p (adj)', dunnett ?
      'Adjusted p-value (Dunnett\'s joint reference distribution). Significant when below α.' :
      'Adjusted p-value (Holm step-down, 1979). Significant when below α.'],
    ['Hedges\' g', 'Hedges\' g effect size. Standardized mean difference, corrected for small-sample bias.'],
  ]);
}

/** CI column tooltip with the per-comparison (not simultaneous) caveat. */
function ciTooltip(alpha: number): string {
  const pct = Math.round((1 - alpha) * 100);
  return `Bound of the ${pct}% confidence interval for the mean difference ` +
    `(per-comparison, not simultaneous).`;
}

// ── Results grid ─────────────────────────────────────────────────

/** One row per comparison (k−1 rows): Conclusion, Group, n, Mean diff, CI, t, df, p, Hedges' g. */
function getControlComparisonsGrid(report: ControlComparisonsReport, factor: DG.Column,
  featureName: string): DG.Grid {
  const comps = report.comparisons;
  const controlLabel = labelOf(factor, report.controlCode);
  const ciLevel = Math.round((1 - report.alpha) * 100);
  const ciCaption = `${ciLevel}% CI`;

  const rawPs = comps.map((c) => c.pValueRaw);
  const adjPs = comps.map((c) => c.pValueAdj);

  const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
    conclusionColumnPerRow(comps.map((c) => c.significant)),
    DG.Column.fromStrings('Group', comps.map((c) => labelOf(factor, c.groupCode))),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', comps.map((c) => c.n)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Mean diff', comps.map((c) => c.meanDiff)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${ciCaption} low`, comps.map((c) => c.ciLow)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${ciCaption} high`, comps.map((c) => c.ciHigh)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', comps.map((c) => c.statistic)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'df', comps.map((c) => c.df)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p (raw)', rawPs),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p (adj)', adjPs),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Hedges\' g', comps.map((c) => c.hedgesG)),
  ]));

  grid.dataFrame.col('p (raw)')!.meta.format = pColumnFormat(rawPs);
  grid.dataFrame.col('p (adj)')!.meta.format = pColumnFormat(adjPs);

  styleConclusionColumn(grid);

  const tooltips = headerTooltips(report, controlLabel, featureName, ciCaption);

  grid.onCellTooltip((cell, x, y) => {
    if (cell.isColHeader) {
      const t = tooltips.get(cell.tableColumn!.name);
      if (t != null) {
        ui.tooltip.show(ui.divV([ui.p(t)]), x, y);
        return true;
      }
    } else if (cell.isTableCell && cell.tableColumn?.name === CONCLUSION_COL_NAME && cell.cell.rowIndex >= 0) {
      const c = comps[cell.cell.rowIndex];
      ui.tooltip.show(conclusionTooltip(c.significant, labelOf(factor, c.groupCode), controlLabel, featureName), x, y);
      return true;
    }
    return false;
  });

  grid.helpUrl = HELP_URL;
  return grid;
} // getControlComparisonsGrid

// ── Layout ───────────────────────────────────────────────────────

/** Box plot docked right + (optionally) the titled results grid docked below. */
function addVisualization(df: DG.DataFrame, factor: DG.Column, feature: DG.Column,
  report: ControlComparisonsReport, showReport: boolean): void {
  const view = grok.shell.getTableView(df.name);
  grok.shell.v = view;

  const controlLabel = labelOf(factor, report.controlCode);

  const chart = DG.Viewer.boxPlot(df, {
    categoryColumnNames: [factor.name],
    valueColumnName: feature.name,
    // A single p-value is meaningless across k ≥ 3 groups; per-group quartiles are useful.
    showPValue: false,
    showStatistics: false,
    description: buildDescription(report, controlLabel),
    descriptionVisibilityMode: 'Always',
    descriptionPosition: 'Top',
    showColorSelector: false,
    showSizeSelector: false,
    autoLayout: false,
  });

  const node = view.dockManager.dock(chart, DG.DOCK_TYPE.RIGHT, null, 'Control comparisons');

  if (!showReport)
    return;

  const grid = getControlComparisonsGrid(report, factor, feature.name);

  // Register the results table in the workspace so users can export / save it.
  grid.dataFrame.name = 'Control comparisons result';
  grok.shell.addTable(grid.dataFrame);

  const title = `Control comparisons: "${feature.name}" vs "${controlLabel}" (${report.method})`;
  view.dockManager.dock(grid, DG.DOCK_TYPE.DOWN, node, title, 0.3);
} // addVisualization

// ── Dialog ───────────────────────────────────────────────────────

function getWarning(msg: string): HTMLElement {
  return ui.divV([
    ui.markdown(`Control comparisons cannot be performed:

    ${msg}`),
    ui.link('Learn more',
      () => window.open('https://en.wikipedia.org/wiki/Dunnett%27s_test', '_blank'),
      'Click to open in a new tab.'),
  ]);
}

/** Run the control-comparisons dialog. */
export function runControlComparisons(): void {
  const df: DG.DataFrame | null = grok.shell.t;

  if (df === null) {
    grok.shell.warning('No dataframe is opened');
    return;
  }

  const columns = df.columns;
  const factorColNames: string[] = [];
  const featureColNames: string[] = [];

  for (const col of columns) {
    if (FEATURE_TYPES.includes(col.type))
      featureColNames.push(col.name);
    else if (FACTOR_TYPES.includes(col.type))
      factorColNames.push(col.name);
  }

  // A category column is usable here if it has at least 2 groups and is not an ID-like (all-unique)
  // column. Used for the applicability check and the default selection; ID columns stay selectable
  // in the dropdown but get flagged (mirrors the ANOVA dialog).
  const validFactorColNames = factorColNames.filter((name) => {
    const col = columns.byName(name);
    return groupCount(col) >= 2 && !isFactorAllUnique(col);
  });

  if (validFactorColNames.length < 1) {
    grok.shell.warning(ui.markdown(`Control comparisons are not applicable to this table:

    - no categorical column with at least two repeated groups`));
    return;
  }

  if (featureColNames.length < 1) {
    grok.shell.warning(ui.markdown(`Control comparisons are not applicable to this table:

    - no numeric feature column (type ${FEATURE_TYPES.join(', ')})`));
    return;
  }

  // Prefer a column with 3+ groups as the default — control comparisons are most meaningful then
  // (2 groups just reduce to a t-test). Fall back to the first valid column otherwise.
  const defaultFactorName = validFactorColNames.find((name) => groupCount(columns.byName(name)) >= 3) ??
    validFactorColNames[0];
  let factor = columns.byName(defaultFactorName);
  let feature = columns.byName(featureColNames[0]);
  let currentMethod: ControlComparisonsMethod = 'Dunnett';
  let significance: number = SIGNIFICANCE.DEFAULT;
  let control: string = groupLabels(factor)[0];

  // --- Method ---
  const methodSource = {method: currentMethod};
  const methodProp = DG.Property.fromOptions({
    name: 'method', caption: 'Method', inputType: 'Radio', choices: METHODS, defaultValue: 'Dunnett',
  });
  const methodInput = ui.input.forProperty(methodProp, methodSource);
  methodInput.onChanged.subscribe(() => {
    currentMethod = (methodSource.method as ControlComparisonsMethod) ?? 'Dunnett';
    updateRunButtonState();
  });
  ui.tooltip.bind(methodInput.captionLabel, () => METHOD_TOOLTIP);

  // --- Category ---
  const factorInput = ui.input.column('Category', {
    table: df,
    value: factor,
    tooltipText: 'Categorical column defining the groups',
    onValueChanged: (col) => {factor = col; rebuildControlChoices(); updateRunButtonState();},
    filter: (col: DG.Column) => factorColNames.includes(col.name),
    nullable: false,
  });

  // Reject ID-like (all-unique) columns via input validation, as in the ANOVA dialog.
  factorInput.addValidator(() => {
    const col = factorInput.value;
    return col != null && isFactorAllUnique(col) ? ALL_UNIQUE_MSG : null;
  });

  // --- Control ---
  const controlInput = ui.input.choice('Control', {
    value: control,
    items: groupLabels(factor),
    tooltipText: 'The reference group every other group is compared against',
    nullable: false,
    onValueChanged: (value) => {control = value as string; updateRunButtonState();},
  });

  function rebuildControlChoices(): void {
    const labels = groupLabels(factor);
    (controlInput as any).items = labels;
    if (!labels.includes(control))
      control = labels[0];
    controlInput.value = control;
  }

  // --- Feature ---
  const featureInput = ui.input.column('Feature', {
    table: df,
    value: feature,
    tooltipText: 'Numeric column to compare across the groups',
    onValueChanged: (col) => {feature = col; updateRunButtonState();},
    filter: (col: DG.Column) => featureColNames.includes(col.name),
    nullable: false,
  });

  // --- Alpha ---
  const signInput = ui.input.float('Alpha', {
    min: SIGNIFICANCE.MIN,
    max: SIGNIFICANCE.MAX,
    value: significance,
    nullable: false,
    tooltipText: 'Significance level',
    onValueChanged: (value) => {significance = value; updateRunButtonState();},
  });

  const fullReportInput = ui.input.bool('Full report', {
    value: true,
    tooltipText: 'Add a table with the full test statistics and per-group conclusions',
  });


  const dlg = ui.dialog({title: 'Control comparisons', helpUrl: HELP_URL});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);

  dlg.addButton('Run', () => {
    dlg.close();
    try {
      const report = controlComparisons(factor!, feature! as NumCol, codeOf(factor!, control),
        binCount(factor!), {method: currentMethod, alpha: significance});
      addVisualization(df, factor!, feature!, report, fullReportInput.value!);
    } catch (error) {
      if (error instanceof Error) {
        grok.shell.warning(getWarning(error.message));
        view.addViewer(DG.VIEWER.BOX_PLOT, {
          categoryColumnNames: [factor!.name],
          valueColumnName: feature!.name,
          showStatistics: true,
        });
      } else
        grok.shell.error('Control comparisons fail: the platform issue');
    }
  }, undefined, 'Perform control comparisons');

  const runBtn = dlg.getButton('Run');

  function updateRunButtonState(): void {
    if (significance <= SIGNIFICANCE.INFIMUM || significance >= SIGNIFICANCE.SUPREMUM) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Alpha must be strictly between 0 and 1');
      return;
    }
    if (factor == null || feature == null) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Select a category column and a feature column.');
      return;
    }

    const nGroups = groupCount(factor);
    if (nGroups < 2) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'The category column needs at least two groups.');
      return;
    }
    // ID-like (all-unique) column: each group has size 1 — block Run (the validator shows why),
    // exactly as in the ANOVA dialog.
    if (isFactorAllUnique(factor)) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, ALL_UNIQUE_MSG);
      return;
    }
    if (control == null || !groupLabels(factor).includes(control)) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Select a control group present in the data.');
      return;
    }

    runBtn.disabled = false;
    ui.tooltip.bind(runBtn, 'Perform control comparisons');
  }

  dlg.add(factorInput)
    .add(controlInput)
    .add(featureInput)
    .add(signInput)
    .add(methodInput)
    .add(fullReportInput);

  updateRunButtonState();
  factorInput.validate();
  dlg.show();

  setTimeout(() => (document.activeElement as HTMLElement | null)?.blur(), 0);
} // runControlComparisons
