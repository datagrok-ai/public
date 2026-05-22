/* eslint-disable max-len */
// Analysis of Variances (ANOVA) - UI

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FactorizedData, oneWayAnova, OneWayAnovaReport, WelchAnova} from './anova-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const ANOVA_HELP_URL = '/help/explore/anova';

/** Significance const */
enum SIGNIFICANCE {
  DEFAULT = 0.05,
  MIN = 0.0001,
  MAX = 0.99,
  INFIMUM = 0,
  SUPREMUM = 1,
};

/** Default names */
enum DEFAULT {
  FACTOR = 'race',
  FEATURE = 'age',
};

type Method = 'Welch' | 'Fisher';

const LEARN_MORE_URL = {
  Fisher: 'https://en.wikipedia.org/wiki/F-test',
  Welch: 'https://en.wikipedia.org/wiki/Welch%27s_t-test',
} as const;


/** Add one-way ANOVA results */
function addVizualization(df: DG.DataFrame, factorsName: string, featuresName: string,
  report: OneWayAnovaReport): void {
  const view = grok.shell.getTableView(df.name);
  grok.shell.v = view;

  const test = report.anovaTable.fStat > report.fCritical;
  const shortConclusion = test ?
    `"${factorsName}" affects the "${featuresName}"` :
    `"${factorsName}" doesn't affect the "${featuresName}"`;

  const chart = DG.Viewer.boxPlot(df, {
    categoryColumnNames: [factorsName],
    valueColumnName: featuresName,
    showPValue: false,
    showStatistics: false,
    description: shortConclusion,
    showColorSelector: false,
    autoLayout: false,
  });

  const node = view.dockManager.dock(chart, DG.DOCK_TYPE.RIGHT, null, 'ANOVA');

  const nullHypoMd = ui.markdown(`**Null Hypothesis:** all group means are equal.`);
  ui.tooltip.bind(nullHypoMd, `The "${factorsName}" factor does not produce a significant difference in the "${featuresName}" feature.`);

  const altHypoMd = ui.markdown(`**Alternative Hypothesis:** at least one group mean differs significantly.`);
  ui.tooltip.bind(altHypoMd, `The "${factorsName}" factor produces a significant difference in the "${featuresName}" feature.`);

  const conclusionMd = ui.markdown(`**Conclusion:** ${test ?
    'significant differences exist between groups.' :
    'no significant differences detected.'}`,
  );

  const tooltipDiv = test ?
    ui.divV([
      ui.p(`Reject the null hypothesis, since F > F-critical:
      ${report.anovaTable.fStat.toFixed(2)} > ${report.fCritical.toFixed(2)}.`),
      ui.h2('There is a significant difference among sample averages.'),
    ]) :
    ui.divV([
      ui.p(`Fail to reject the null hypothesis, since F < F-critical:
      ${report.anovaTable.fStat.toFixed(2)} < ${report.fCritical.toFixed(2)}.`),
      ui.h2('There is no significant difference among sample averages.'),
    ]);

  ui.tooltip.bind(conclusionMd, () => tooltipDiv);

  const divResult = ui.divV([
    nullHypoMd,
    altHypoMd,
    conclusionMd,
    ui.link('Learn more',
      () => window.open(LEARN_MORE_URL[report.method], '_blank'),
      'Click to open in a new tab.',
    ),
  ]);

  const analysisTitle = report.method === 'Welch' ?
    "One-Way ANOVA (Welch's)" :
    "One-Way ANOVA (Fisher's)";
  const reportViewer = getAnovaGrid(report);
  const tabControl = ui.tabControl({
    'Analysis': ui.panel([ui.h2(analysisTitle), reportViewer.root]),
    'Conclusion': ui.panel([divResult]),
  });

  ui.tooltip.bind(tabControl.getPane('Analysis').header, 'ANOVA results summary.');
  ui.tooltip.bind(tabControl.getPane('Conclusion').header, 'Null hypothesis testing.');


  view.dockManager.dock(tabControl.root, DG.DOCK_TYPE.DOWN, node, '', 0.25);

  reportViewer.root.style.width = '100%';
} // addVizualization

/** Create grid with one-way ANOVA results, dispatched by method. */
function getAnovaGrid(report: OneWayAnovaReport): DG.Grid {
  if (report.method === 'Fisher')
    return getFisherGrid(report);
  return getWelchGrid(report.anovaTable, report.fCritical, report.significance);
}

/** Classical Fisher ANOVA table (3x7: Between/Within/Total x SS/DF/MS/F/F-critical/p-value). */
function getFisherGrid(report: OneWayAnovaReport & {method: 'Fisher'}): DG.Grid {
  const anova = report.anovaTable;

  const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Source of variance', ['Between groups', 'Within groups', 'Total']),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'SS', [anova.ssBn, anova.ssWn, anova.ssTot]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'DF', [anova.dfBn, anova.dfWn, anova.dfTot]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'MS', [anova.msBn, anova.msWn, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F', [anova.fStat, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F-critical', [report.fCritical, null, null]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p-value', [anova.pValue, null, null]),
  ]));

  const tooltip = new Map([
    ['Source of variance', 'List of the explored variation sources.'],
    ['SS', 'Sum of squares (SS). Measure of total variation in the data.'],
    ['DF', 'Degrees of freedom (DF). Number of independent values that can vary.'],
    ['MS', 'Mean square (MS). Sum of squares divided by degrees of freedom.'],
    ['F', 'F-statistics (F). Ratio of between-group to within-group variance.'],
    ['F-critical', `${report.significance}-critical value of F-statistics.`],
    ['p-value', `Probability of observing this result if groups have equal means.`],
  ]);

  grid.onCellTooltip(function(cell, x, y) {
    if (cell.isColHeader) {
      ui.tooltip.show(ui.divV([ui.p(tooltip.get(cell.tableColumn!.name)!)]), x, y);
      return true;
    }
  });

  grid.helpUrl = ANOVA_HELP_URL;

  return grid;
} // getFisherGrid

/** Welch ANOVA results (1x6: Source / F / df₁ / df₂ / F-critical / p-value).
 *  Welch's W-test does not have a SS/MS decomposition. p-value is rendered
 *  in APA style ('< .001', '.04') via a STRING column. */
function getWelchGrid(anova: WelchAnova, fCritical: number, significance: number): DG.Grid {
  const DF1 = 'df₁';
  const DF2 = 'df₂';

  const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Source of variance', ['Between groups']),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F', [anova.fStat]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, DF1, [anova.dfBn]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, DF2, [anova.dfWn]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'F-critical', [fCritical]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p-value', [anova.pValue]),
  ]));

  const tooltip = new Map([
    ['Source of variance', 'List of the explored variation sources.'],
    ['F', `F-statistic. Ratio of weighted between-group variance to within-group variance (Welch's W).`],
    [DF1, 'Numerator degrees of freedom (k − 1, where k is the number of groups).'],
    [DF2, 'Welch–Satterthwaite-adjusted denominator degrees of freedom. Fractional by design.'],
    ['F-critical', `${significance}-critical value of F-statistic with Welch's df.`],
    ['p-value', 'Probability of observing this result if group means are equal.'],
  ]);

  grid.onCellTooltip(function(cell, x, y) {
    if (cell.isColHeader) {
      ui.tooltip.show(ui.divV([ui.p(tooltip.get(cell.tableColumn!.name)!)]), x, y);
      return true;
    }
  });

  grid.helpUrl = ANOVA_HELP_URL;

  return grid;
} // getWelchGrid

/** Return warning div */
function getWarning(msg: string): HTMLElement {
  return ui.divV([
    ui.markdown(`ANOVA cannot be performed:

    ${msg}`),
    ui.link('Learn more',
      () => window.open('https://en.wikipedia.org/wiki/Analysis_of_variance#Assumptions', '_blank'),
      'Click to open in a new tab',
    ),
  ]);
}

/** Run one-way analysis of variances */
export function runOneWayAnova(): void {
  /** current dataframe */
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

  const factorColsCount = factorColNames.length;
  if (factorColsCount < 1) {
    grok.shell.warning(ui.markdown(`No acceptable factor columns:

    - type: ${FACTOR_TYPES.join(', ')}
    - at least two categories`,
    ));
    return;
  };

  let factor = df.col(DEFAULT.FACTOR);

  if (factor === null) {
    let minIdx = 0;
    let minCount = columns.byName(factorColNames[0]).categories.length;
    let current: number;

    for (let i = 1; i < factorColsCount; ++i) {
      current = columns.byName(factorColNames[i]).categories.length;
      if (current < minCount) {
        minCount = current;
        minIdx = i;
      }
    }

    factor = columns.byName(factorColNames[minIdx]);
  }

  if (featureColNames.length < 1) {
    grok.shell.warning(ui.markdown(`No acceptable feature columns:

    - type: ${FEATURE_TYPES.join(', ')}`,
    ));
    return;
  }

  const factorInput = ui.input.column('Category', {
    table: df,
    value: factor,
    tooltipText: 'Column with factor values',
    onValueChanged: (col) => {factor = col; updateRunButtonState();},
    filter: (col: DG.Column) => factorColNames.includes(col.name),
    nullable: false,
  });

  let feature = df.col(DEFAULT.FEATURE);
  if (feature === null)
    feature = columns.byName(featureColNames[0]);

  const featureInput = ui.input.column('Feature', {
    table: df,
    value: feature,
    tooltipText: 'Column with feature values',
    onValueChanged: (col) => {feature = col; updateRunButtonState();},
    filter: (col: DG.Column) => featureColNames.includes(col.name),
    nullable: false,
  });

  let currentMethod: Method = 'Welch';
  const methodSource = {method: currentMethod};
  const methodProp = DG.Property.fromOptions({
    name: 'method',
    caption: 'Method',
    inputType: 'Radio',
    choices: ['Welch', 'Fisher'],
    defaultValue: 'Welch',
  });
  const methodInput = ui.input.forProperty(methodProp, methodSource);
  methodInput.onChanged.subscribe(() => {
    currentMethod = (methodSource.method as Method) ?? 'Welch';
    updateRunButtonState();
  });

  const methodTooltip = ui.markdown(
    'Set the method for analysis:\n\n' +
    '* **Welch** — robust to **unequal variances** across groups. ' +
      'Recommended as the default when group variances are unknown or differ.\n\n' +
    '* **Fisher** — classical ANOVA. Assumes **equal variances** across groups; ' +
      'more powerful when this assumption holds.\n\n',
  );
  ui.tooltip.bind(methodInput.captionLabel, () => methodTooltip);

  const FISHER_UNEQUAL_VAR_MSG =
    'Variances differ significantly between groups. ' +
    `Fisher's ANOVA requires equal variances — switch Method to Welch.`;
  const fisherWarningIcon = ui.iconFA('info-circle', null, FISHER_UNEQUAL_VAR_MSG);
  fisherWarningIcon.style.color = 'var(--red-3, #EB6767)';
  fisherWarningIcon.style.marginLeft = '12px';
  fisherWarningIcon.style.display = 'none';
  methodInput.root.append(fisherWarningIcon);

  let significance = SIGNIFICANCE.DEFAULT;
  const signInput = ui.input.float('Alpha', {
    min: SIGNIFICANCE.MIN,
    max: SIGNIFICANCE.MAX,
    value: significance,
    nullable: false,
    tooltipText: 'Significance level',
    onValueChanged: (value) => {
      significance = value;
      updateRunButtonState();
    },
  });

  const dlg = ui.dialog({title: 'ANOVA', helpUrl: ANOVA_HELP_URL});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);
  dlg.addButton('Run', () => {
    dlg.close();

    try {
      const res = oneWayAnova(factor!, feature!, significance, {
        method: currentMethod,
        toValidate: false,
      });
      addVizualization(df, factor!.name, feature!.name, res);
    } catch (error) {
      if (error instanceof Error) {
        grok.shell.warning(getWarning(error.message));

        view.addViewer(DG.VIEWER.BOX_PLOT, {
          categoryColumnNames: [factor!.name],
          valueColumnName: feature!.name,
          showStatistics: false,
          showPValue: false,
        });
      } else
        grok.shell.error('ANOVA fails: the platform issue');
    }
  }, undefined, 'Perform analysis of variances');

  const runBtn = dlg.getButton('Run');

  function updateRunButtonState(): void {
    fisherWarningIcon.style.display = 'none';

    if (significance <= SIGNIFICANCE.INFIMUM || significance >= SIGNIFICANCE.SUPREMUM) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, 'Alpha must be strictly between 0 and 1.');
      return;
    }

    let varEqual: boolean;
    try {
      const uniqueCount = factor!.stats.uniqueCount;
      if (uniqueCount < 2)
        throw new Error('At least two categories required.');
      const factorized = new FactorizedData(factor!, feature!, uniqueCount);
      varEqual = factorized.areVarsEqual(significance);
    } catch (err) {
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, (err as Error).message);
      return;
    }

    if (currentMethod === 'Fisher' && !varEqual) {
      fisherWarningIcon.style.display = '';
      runBtn.disabled = true;
      ui.tooltip.bind(runBtn, FISHER_UNEQUAL_VAR_MSG);
      return;
    }

    runBtn.disabled = false;
    ui.tooltip.bind(runBtn, 'Perform analysis of variances');
  }

  dlg.add(factorInput)
    .add(featureInput)
    .add(methodInput)
    .add(signInput);

  updateRunButtonState();

  dlg.show();

  // Strip auto-focus from the Method radio: when the dialog opens, focus
  // lands on the first focusable input and a focus-ring is drawn around
  // the selected Welch option — visually distracting since the user
  // hasn't interacted yet.
  setTimeout(() => (document.activeElement as HTMLElement | null)?.blur(), 0);
} // runOneWayAnova
