/* eslint-disable max-len */
// Two-sample t-test - UI.
//
// Mirrors the ANOVA dialog (`runOneWayAnova`). Reuses the ANOVA factorization and
// stat primitives, but ships its own results renderer (the unification of the ANOVA
// and t-test renderers into a single two-level discriminated union is a follow-up).

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {factorize, groupStats, areVarsEqual, GroupStats, CatCol, NumCol} from '../anova/anova-tools';
import {twoSampleTTest, TwoSampleTTest, TTestMethod} from './ttest-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const T_TEST_HELP_URL = '/help/explore/t-test';

/** Significance const. */
enum SIGNIFICANCE {
  DEFAULT = 0.05,
  MIN = 0.0001,
  MAX = 0.99,
  INFIMUM = 0,
  SUPREMUM = 1,
};

const LEARN_MORE_URL = {
  Student: 'https://en.wikipedia.org/wiki/Student%27s_t-test',
  Welch: 'https://en.wikipedia.org/wiki/Welch%27s_t-test',
} as const;

const STUDENT_UNEQUAL_VAR_MSG =
  'Variances differ significantly between the two groups. ' +
  `Student's t-test assumes equal variances — switch Method to Welch.`;

/** Number of factorization bins for a categorical column (upper bound on category codes). */
function binCount(factor: DG.Column): number {
  return factor.type === DG.COLUMN_TYPE.BOOL ? 2 : factor.categories.length;
}

/** Number of non-empty groups in a categorical column (distinct non-null values). */
function groupCount(factor: DG.Column): number {
  return factor.stats.uniqueCount;
}

/** Non-empty groups (with their category codes) produced by factorizing `feature` by `factor`. */
function presentGroups(factor: CatCol, feature: NumCol): {code: number, stats: GroupStats}[] {
  const f = factorize(factor, feature, binCount(factor));
  const present: {code: number, stats: GroupStats}[] = [];
  for (let i = 0; i < f.catCount; ++i) {
    if (f.sizes[i] > 0)
      present.push({code: i, stats: groupStats(f, i)});
  }
  return present;
}

/** Human-readable label for a category code. */
function labelOf(factor: DG.Column, code: number): string {
  if (factor.type === DG.COLUMN_TYPE.BOOL)
    return code === 1 ? 'true' : 'false';
  return String(factor.categories[code]);
}

/** Format a number for display (compact scientific notation for tiny magnitudes). */
function fmt(x: number): string {
  if (!Number.isFinite(x))
    return String(x);
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-4 || a >= 1e6))
    return x.toExponential(2);
  return x.toFixed(4);
}

/** Qualitative magnitude of a standardized effect size (Cohen's conventions). */
function effectMagnitude(g: number): string {
  const a = Math.abs(g);
  if (a < 0.2) return 'negligible';
  if (a < 0.5) return 'small';
  if (a < 0.8) return 'medium';
  return 'large';
}

/** One-row results grid: t, df, p-value, mean difference, 95% CI, effect size. */
function getTTestGrid(res: TwoSampleTTest, factorName: string): DG.Grid {
  const ciCaption = `${Math.round((1 - res.alpha) * 100)}% CI`;
  const effectCaption = res.method === 'Welch' ? 'g_s*' : 'Hedges\' g';

  const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 't', [res.t]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'df', [res.df]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p-value', [res.pValue]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Mean difference', [res.meanDiff]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${ciCaption} low`, [res.ciLow]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${ciCaption} high`, [res.ciHigh]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Cohen\'s d', [res.cohenD]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, effectCaption, [res.hedgesG]),
  ]));

  const tooltip = new Map<string, string>([
    ['t', 'Two-sample t-statistic (signed; group 1 − group 0).'],
    ['df', res.method === 'Welch' ?
      'Welch–Satterthwaite degrees of freedom (fractional by design).' :
      'Degrees of freedom (n0 + n1 − 2).'],
    ['p-value', 'Two-sided probability of observing this result if the group means are equal.'],
    ['Mean difference', `Difference of the two group means (group 1 − group 0) in the "${factorName}" split.`],
    [`${ciCaption} low`, `Lower bound of the ${ciCaption} for the mean difference.`],
    [`${ciCaption} high`, `Upper bound of the ${ciCaption} for the mean difference.`],
    ['Cohen\'s d', res.method === 'Welch' ?
      'Standardized mean difference (d_s) using the non-pooled SD √((v0+v1)/2).' :
      'Standardized mean difference using the pooled SD.'],
    [effectCaption, res.method === 'Welch' ?
      'Bias-corrected effect size (Hedges\' g_s*), non-pooled SD.' :
      'Bias-corrected effect size (Hedges\' g), pooled SD.'],
  ]);

  grid.onCellTooltip(function(cell, x, y) {
    if (cell.isColHeader) {
      const t = tooltip.get(cell.tableColumn!.name);
      if (t != null) {
        ui.tooltip.show(ui.divV([ui.p(t)]), x, y);
        return true;
      }
    }
    return false;
  });

  grid.helpUrl = T_TEST_HELP_URL;
  return grid;
} // getTTestGrid

/** Lay out the t-test results: box plot + Analysis/Conclusion tab control. */
function addVisualization(df: DG.DataFrame, factor: DG.Column, feature: DG.Column,
  res: TwoSampleTTest): void {
  const view = grok.shell.getTableView(df.name);
  grok.shell.v = view;

  const significant = res.pValue < res.alpha;
  const label0 = labelOf(factor, res.code0);
  const label1 = labelOf(factor, res.code1);

  const shortConclusion = significant ?
    `"${feature.name}" differs significantly between "${factor.name}" groups` :
    `"${feature.name}" shows no significant difference between "${factor.name}" groups`;

  const chart = DG.Viewer.boxPlot(df, {
    categoryColumnNames: [factor.name],
    valueColumnName: feature.name,
    showPValue: true,
    showStatistics: true,
    description: shortConclusion,
    showColorSelector: false,
    autoLayout: false,
  });

  const node = view.dockManager.dock(chart, DG.DOCK_TYPE.RIGHT, null, 'T-Test');

  // Conclusion: "p = …, effect = …, 95% CI […]" — not a bare reject/fail.
  const effectCaption = res.method === 'Welch' ? 'g_s*' : 'Hedges\' g';
  const ciCaption = `${Math.round((1 - res.alpha) * 100)}% CI`;

  const verdictMd = ui.markdown(`**Result:** ${significant ?
    `a significant difference between "${label0}" and "${label1}".` :
    `no significant difference between "${label0}" and "${label1}".`}`);

  const summaryMd = ui.markdown(
    `**p** = ${fmt(res.pValue)} ${significant ? '<' : '≥'} α = ${res.alpha}  ·  ` +
    `**effect** (${effectCaption}) = ${fmt(res.hedgesG)} (${effectMagnitude(res.hedgesG)})  ·  ` +
    `**${ciCaption}** for the mean difference = [${fmt(res.ciLow)}, ${fmt(res.ciHigh)}]`,
  );

  const meanDiffMd = ui.markdown(
    `Mean difference (mean("${label1}") − mean("${label0}")) = ${fmt(res.meanDiff)} ` +
    `(t = ${fmt(res.t)}, df = ${fmt(res.df)}).`,
  );

  const divResult = ui.divV([
    verdictMd,
    summaryMd,
    meanDiffMd,
    ui.link('Learn more',
      () => window.open(LEARN_MORE_URL[res.method], '_blank'),
      'Click to open in a new tab.',
    ),
  ]);

  const analysisTitle = res.method === 'Welch' ?
    'Two-Sample t-test (Welch\'s)' :
    'Two-Sample t-test (Student\'s)';
  const reportViewer = getTTestGrid(res, factor.name);

  const tabControl = ui.tabControl({
    'Analysis': ui.panel([ui.h2(analysisTitle), reportViewer.root]),
    'Conclusion': ui.panel([divResult]),
  });

  ui.tooltip.bind(tabControl.getPane('Analysis').header, 't-test results summary.');
  ui.tooltip.bind(tabControl.getPane('Conclusion').header, 'Significance and effect size.');

  view.dockManager.dock(tabControl.root, DG.DOCK_TYPE.DOWN, node, '', 0.25);

  reportViewer.root.style.width = '100%';
} // addVisualization

/** Warning div shown when the t-test cannot be performed. */
function getWarning(msg: string): HTMLElement {
  return ui.divV([
    ui.markdown(`The t-test cannot be performed:

    ${msg}`),
    ui.link('Learn more',
      () => window.open(LEARN_MORE_URL.Welch, '_blank'),
      'Click to open in a new tab.',
    ),
  ]);
}

/** Run the two-sample t-test dialog. */
export function runTwoSampleTTest(): void {
  const df: DG.DataFrame | null = grok.shell.t;

  if (df === null) {
    grok.shell.warning('No dataframe is opened');
    return;
  }

  const columns = df.columns;
  const factorColNames = [] as string[];
  const featureColNames = [] as string[];

  for (const col of columns) {
    if (FEATURE_TYPES.includes(col.type))
      featureColNames.push(col.name);
    else if (FACTOR_TYPES.includes(col.type))
      factorColNames.push(col.name);
  }

  if (factorColNames.length < 1) {
    grok.shell.warning(ui.markdown(`The t-test is not applicable to this table:

    - no categorical column (type ${FACTOR_TYPES.join(', ')}) to define groups`,
    ));
    return;
  }

  if (featureColNames.length < 1) {
    grok.shell.warning(ui.markdown(`The t-test is not applicable to this table:

    - no numeric feature column (type ${FEATURE_TYPES.join(', ')})`,
    ));
    return;
  }

  // A categorical column is valid for the two-sample t-test only if it has exactly two groups.
  const validFactorColNames = factorColNames.filter((name) => groupCount(columns.byName(name)) === 2);

  if (validFactorColNames.length < 1) {
    grok.shell.warning(ui.markdown(`The t-test is not applicable to this table:

    - no categorical column with exactly two groups
    - the two-sample t-test compares exactly two groups — use ANOVA for three or more`,
    ));
    return;
  }

  // Default to the first valid (two-group) category column.
  let factor = columns.byName(validFactorColNames[0]);

  let feature = columns.byName(featureColNames[0]);

  let currentMethod: TTestMethod = 'Welch';
  let significance: number = SIGNIFICANCE.DEFAULT;

  // --- Method (Welch / Student) ---
  const methodSource = {method: currentMethod};
  const methodProp = DG.Property.fromOptions({
    name: 'method',
    caption: 'Method',
    inputType: 'Radio',
    choices: ['Welch', 'Student'],
    defaultValue: 'Welch',
  });
  const methodInput = ui.input.forProperty(methodProp, methodSource);
  methodInput.onChanged.subscribe(() => {
    currentMethod = (methodSource.method as TTestMethod) ?? 'Welch';
    updateRunButtonState();
  });

  const methodTooltip = ui.markdown(
    'Set the method for the test:\n\n' +
    '* **Welch** — robust to **unequal variances**. Recommended default.\n\n' +
    '* **Student** — assumes **equal variances**; more powerful when that holds.\n\n',
  );
  ui.tooltip.bind(methodInput.captionLabel, () => methodTooltip);

  const studentWarningIcon = ui.iconFA('info-circle', null, STUDENT_UNEQUAL_VAR_MSG);
  studentWarningIcon.style.color = 'var(--red-3, #EB6767)';
  studentWarningIcon.style.marginLeft = '12px';
  studentWarningIcon.style.display = 'none';
  methodInput.root.append(studentWarningIcon);

  // --- Category (group) ---
  const factorInput = ui.input.column('Category', {
    table: df,
    value: factor,
    tooltipText: 'Categorical column defining the two groups',
    onValueChanged: (col) => {factor = col; updateRunButtonState();},
    filter: (col: DG.Column) => factorColNames.includes(col.name),
    nullable: false,
  });

  // Red indicator on the Category input for an invalid category (≠ 2 groups, or a
  // degenerate group) — by analogy with the Student/Welch indicator on Method.
  const factorWarningIcon = ui.iconFA('info-circle', null);
  factorWarningIcon.style.color = 'var(--red-3, #EB6767)';
  factorWarningIcon.style.marginLeft = '12px';
  factorWarningIcon.style.display = 'none';
  factorInput.root.append(factorWarningIcon);

  // --- Feature ---
  const featureInput = ui.input.column('Feature', {
    table: df,
    value: feature,
    tooltipText: 'Numeric column to compare across the two groups',
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

  const dlg = ui.dialog({title: 't-test', helpUrl: T_TEST_HELP_URL});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);

  dlg.addButton('Run', () => {
    dlg.close();
    try {
      const res = twoSampleTTest(factor!, feature!, binCount(factor!), {
        method: currentMethod,
        alpha: significance,
      });
      addVisualization(df, factor!, feature!, res);
    } catch (error) {
      if (error instanceof Error) {
        grok.shell.warning(getWarning(error.message));
        view.addViewer(DG.VIEWER.BOX_PLOT, {
          categoryColumnNames: [factor!.name],
          valueColumnName: feature!.name,
          showStatistics: true,
          showPValue: true,
        });
      } else
        grok.shell.error('t-test fails: the platform issue');
    }
  }, undefined, 'Perform two-sample t-test');

  const runBtn = dlg.getButton('Run');

  function updateRunButtonState(): void {
    studentWarningIcon.style.display = 'none';
    factorWarningIcon.style.display = 'none';

    /** Flag an invalid Category: show the red icon + disable Run, both with `msg`. */
    const flagCategory = (msg: string): void => {
      factorWarningIcon.style.display = '';
      ui.tooltip.bind(factorWarningIcon, msg);
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, msg);
    };

    // 1. Alpha range.
    if (significance <= SIGNIFICANCE.INFIMUM || significance >= SIGNIFICANCE.SUPREMUM) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Alpha must be strictly between 0 and 1.');
      return;
    }

    if (factor == null || feature == null) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Select a category column and a feature column.');
      return;
    }

    let present: {code: number, stats: GroupStats}[];
    try {
      present = presentGroups(factor, feature);
    } catch (err) {
      flagCategory((err as Error).message);
      return;
    }

    // 2. Number of groups must be exactly 2.
    if (present.length > 2) {
      flagCategory(`Column has ${present.length} groups (3+). Use ANOVA.`);
      return;
    }
    if (present.length < 2) {
      flagCategory('Exactly 2 groups required — the column has fewer than two.');
      return;
    }

    // 3. Degenerate groups (variance cannot be estimated).
    if (present[0].stats.size < 2 || present[1].stats.size < 2) {
      flagCategory('Each group needs at least 2 observations to estimate variance.');
      return;
    }

    // 4. Student + unequal variances → suggest Welch.
    if (currentMethod === 'Student' && !areVarsEqual(present[0].stats, present[1].stats, significance)) {
      studentWarningIcon.style.display = '';
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, STUDENT_UNEQUAL_VAR_MSG);
      return;
    }

    runBtn.disabled = false;
    ui.tooltip.bind(runBtn, 'Perform two-sample t-test');
  }

  dlg.add(factorInput)
    .add(featureInput)
    .add(signInput)
    .add(methodInput);

  updateRunButtonState();

  dlg.show();

  // Strip auto-focus from the first input to avoid a distracting focus ring.
  setTimeout(() => (document.activeElement as HTMLElement | null)?.blur(), 0);
} // runTwoSampleTTest
