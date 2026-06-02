/* eslint-disable max-len */
// Two-sample t-test - UI.
//
// Mirrors the ANOVA dialog (`runOneWayAnova`). Reuses the ANOVA factorization and
// stat primitives, but ships its own results renderer.

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {factorize, groupStats, areVarsEqual, GroupStats, CatCol, NumCol} from '../anova/anova-tools';
import {twoSampleTTest, TwoSampleTTest, TTestMethod} from './ttest-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const T_TEST_HELP_URL = '/help/explore/group-comparison#t-test';

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

/** Qualitative magnitude of a standardized effect size (Cohen's conventions; domain may dictate different cutoffs). */
function effectMagnitude(g: number): string {
  const a = Math.abs(g);
  if (a < 0.2) return 'negligible';
  if (a < 0.5) return 'small';
  if (a < 0.8) return 'medium';
  return 'large';
}

// --- Box plot description (see boxplot-description-spec.md) ---

/** Separator between the phrase and the p-value in the box plot description (placed directly after the phrase). */
const DESCRIPTION_SEPARATOR = ',';
/** Above this total line length the description switches to the long (category-free) form. */
const DESCRIPTION_LENGTH_THRESHOLD = 70;

/** Format a p-value for the box plot description: clamp at the extremes, else round to 3 digits without trailing zeros. */
function formatPForDescription(p: number): string {
  if (p < 0.001)
    return 'p < 0.001';
  if (p > 0.99)
    return 'p > 0.99';
  // Round to 3 decimals; Number() drops trailing zeros (0.420 → 0.42, 0.005 → 0.005).
  return `p = ${Number(p.toFixed(3))}`;
}

/**
 * Build the box plot description line per boxplot-description-spec.md.
 *
 * Short form ("<value>" differs between "<catA>" and "<catB>") is used when the whole line fits
 * within the length threshold and both category names are non-empty; otherwise the long, category-free
 * form ("<factor>" affects the "<value>") is used.
 */
function buildDescription(factorName: string, featureName: string, label0: string, label1: string,
  significant: boolean, p: number): string {
  const pStr = formatPForDescription(p);

  const shortPhrase = significant ?
    `"${featureName}" differs between "${label0}" and "${label1}"` :
    `"${featureName}" doesn't differ between "${label0}" and "${label1}"`;
  const shortLine = `${shortPhrase}${DESCRIPTION_SEPARATOR} ${pStr}`;

  const labelsPresent = label0.length > 0 && label1.length > 0;
  if (labelsPresent && shortLine.length <= DESCRIPTION_LENGTH_THRESHOLD)
    return shortLine;

  // Long form: note the deliberate article ("affects the") absent from the short form.
  const longPhrase = significant ?
    `"${factorName}" affects the "${featureName}"` :
    `"${factorName}" doesn't affect the "${featureName}"`;
  return `${longPhrase}${DESCRIPTION_SEPARATOR} ${pStr}`;
}

/** One-row results grid: t, df, p-value, mean difference, 95% CI, effect size. */
function getTTestGrid(res: TwoSampleTTest, label0: string, label1: string): DG.Grid {
  const ciLevel = Math.round((1 - res.alpha) * 100);
  const ciCaption = `${ciLevel}% CI`;
  // Single user-facing header for both methods; tooltip explains the pooled vs non-pooled distinction.
  const effectCaption = 'Hedges\' g';
  const diffDescription = `mean("${label1}") − mean("${label0}")`;

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
    ['t', `Two-sample t-statistic. Positive when mean("${label1}") > mean("${label0}").`],
    ['df', res.method === 'Welch' ?
      'Welch–Satterthwaite degrees of freedom. Non-integer values are normal for Welch.' :
      'Degrees of freedom (n0 + n1 − 2).'],
    ['p-value',
      'Probability of seeing a difference at least as extreme as the observed one, ' +
      'assuming the two group means are actually equal. Smaller p = stronger evidence against equality.'],
    ['Mean difference', `${diffDescription}, in the original units of the feature.`],
    [`${ciCaption} low`, `Lower bound of the ${ciCaption} (confidence interval) for ${diffDescription}.`],
    [`${ciCaption} high`, `Upper bound of the ${ciCaption} (confidence interval) for ${diffDescription}.`],
    ['Cohen\'s d', res.method === 'Welch' ?
      'Standardized mean difference. For Welch, uses the average of the two group SDs (non-pooled, per Delacre et al. 2021).' :
      'Standardized mean difference, using the pooled SD across the two groups.'],
    [effectCaption, res.method === 'Welch' ?
      `Hedges' g_s* — Cohen's d corrected for small-sample bias. Companion effect size to Welch's t-test (Delacre et al. 2021). ` +
      `Rules of thumb (Cohen 1988): |g| < 0.2 negligible, < 0.5 small, < 0.8 medium, ≥ 0.8 large.` :
      `Hedges' g — Cohen's d corrected for small-sample bias. ` +
      `Rules of thumb (Cohen 1988): |g| < 0.2 negligible, < 0.5 small, < 0.8 medium, ≥ 0.8 large.`],
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

/** Lay out the t-test results: box plot + Analysis/Conclusion tab control (mirrors ANOVA). */
function addVisualization(df: DG.DataFrame, factor: DG.Column, feature: DG.Column,
  res: TwoSampleTTest): void {
  const view = grok.shell.getTableView(df.name);
  grok.shell.v = view;

  const significant = res.pValue < res.alpha;
  const label0 = labelOf(factor, res.code0);
  const label1 = labelOf(factor, res.code1);

  const description = buildDescription(factor.name, feature.name, label0, label1, significant, res.pValue);

  const chart = DG.Viewer.boxPlot(df, {
    categoryColumnNames: [factor.name],
    valueColumnName: feature.name,
    // No on-plot statistics or p-value (as in ANOVA); all numbers live in the Analysis/Conclusion tabs.
    showPValue: false,
    showStatistics: false,
    description: description,
    // Always show the description, pinned to the top of the box plot.
    descriptionVisibilityMode: 'Always',
    descriptionPosition: 'Top',
    showColorSelector: false,
    showSizeSelector: false,
    autoLayout: false,
  });

  const node = view.dockManager.dock(chart, DG.DOCK_TYPE.RIGHT, null, 'T-Test');

  const effectCaption = res.method === 'Welch' ? 'Hedges\' g_s*' : 'Hedges\' g';
  const ciCaption = `${Math.round((1 - res.alpha) * 100)}% CI`;

  const nullHypoMd = ui.markdown(`**Null Hypothesis:** the two group means are equal.`);
  ui.tooltip.bind(nullHypoMd, `The mean of "${feature.name}" is the same for "${label0}" and "${label1}".`);

  const altHypoMd = ui.markdown(`**Alternative Hypothesis:** the two group means differ.`);
  ui.tooltip.bind(altHypoMd, `The mean of "${feature.name}" differs between "${label0}" and "${label1}".`);

  const conclusionMd = ui.markdown(`**Conclusion:** ${significant ?
    'a significant difference between the group means.' :
    'no significant difference between the group means.'}`,
  );

  const tooltipDiv = significant ?
    ui.divV([
      ui.p(`Reject the null hypothesis, since p < α: ${fmt(res.pValue)} < ${res.alpha}.`),
      ui.p(`Effect size (${effectCaption}) = ${fmt(res.hedgesG)} (${effectMagnitude(res.hedgesG)}); ` +
        `${ciCaption} for the mean difference = [${fmt(res.ciLow)}, ${fmt(res.ciHigh)}].`),
      ui.h2('There is a significant difference between the group means.'),
    ]) :
    ui.divV([
      ui.p(`Fail to reject the null hypothesis, since p ≥ α: ${fmt(res.pValue)} ≥ ${res.alpha}.`),
      ui.p(`Effect size (${effectCaption}) = ${fmt(res.hedgesG)} (${effectMagnitude(res.hedgesG)}); ` +
        `${ciCaption} for the mean difference = [${fmt(res.ciLow)}, ${fmt(res.ciHigh)}].`),
      ui.h2('There is no significant difference between the group means.'),
    ]);

  ui.tooltip.bind(conclusionMd, () => tooltipDiv);

  const divResult = ui.divV([
    nullHypoMd,
    altHypoMd,
    conclusionMd,
    ui.link('Learn more',
      () => window.open(LEARN_MORE_URL[res.method], '_blank'),
      'Click to open in a new tab.',
    ),
  ]);

  const analysisTitle = res.method === 'Welch' ?
    'Two-Sample t-test (Welch\'s)' :
    'Two-Sample t-test (Student\'s)';
  const reportViewer = getTTestGrid(res, label0, label1);

  const tabControl = ui.tabControl({
    'Analysis': ui.panel([ui.h2(analysisTitle), reportViewer.root]),
    'Conclusion': ui.panel([divResult]),
  });

  ui.tooltip.bind(tabControl.getPane('Analysis').header, 't-test results summary.');
  ui.tooltip.bind(tabControl.getPane('Conclusion').header, 'Null hypothesis testing.');

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
    '* **Welch** — does **not** assume equal group variances. Recommended default; ' +
    'controls error rates correctly whether variances are equal or not, ' +
    'with negligible power loss when they are.\n\n' +
    '* **Student** — assumes **equal variances** across the two groups. ' +
    'Slightly more powerful when that assumption truly holds; otherwise its p-value is unreliable.\n\n',
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

  const dlg = ui.dialog({title: 'Two-sample t-test', helpUrl: T_TEST_HELP_URL});
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
