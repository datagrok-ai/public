// Analysis of Variances (ANOVA): UI

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {oneWayAnova, OneWayAnovaReport} from './anova-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const ANOVA_HELP_URL = '/help/explore/anova';

/** Significance const */
enum SIGNIFICANCE {
  DEFAULT = 0.05,
  MIN = 0.01,
  MAX = 0.99,
  INFIMUM = 0,
  SUPREMUM = 1,
};

/** Default names */
enum DEFAULT {
  FACTOR = 'race',
  FEATURE = 'age',
};

/** Add one-way ANOVA results */
function addVizualization(df: DG.DataFrame, factorsName: string, featuresName: string, report: OneWayAnovaReport) {
  const test = report.anovaTable.fStat > report.fCritical;

  const shortConclusion = test ?
    `"${factorsName}" affects the "${featuresName}"` :
    `"${factorsName}" doesn't affect the "${featuresName}"`;

  const view = grok.shell.getTableView(df.name);
  const boxPlot = DG.Viewer.boxPlot(df, {
    categoryColumnNames: [factorsName],
    valueColumnName: featuresName,
    showPValue: false,
    showStatistics: false,
    description: shortConclusion,
    showColorSelector: false,
  });
  const boxPlotNode = view.dockManager.dock(boxPlot.root, DG.DOCK_TYPE.RIGHT, null, 'ANOVA');

  const hypoMd = ui.markdown(`**H0:** the "${factorsName}" 
  factor does not produce a significant difference in the "${featuresName}" feature.`);
  ui.tooltip.bind(hypoMd, 'Null hypothesis');

  const testMd = ui.markdown(`**Test result:** ${test ?
    'means differ significantly.' :
    'means do not differ significantly.'}`,
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

  ui.tooltip.bind(testMd, () => tooltipDiv);

  const divResult = ui.divV([
    hypoMd,
    testMd,
    ui.link('Learn more',
      () => window.open('https://en.wikipedia.org/wiki/F-test', '_blank'),
      'Click to open in a new tab',
    ),
  ]);
  divResult.style.marginLeft = '20px';

  const hypoNode = grok.shell.dockManager.dock(divResult, DG.DOCK_TYPE.DOWN, boxPlotNode, 'F-test', 0.3);

  const reportViewer = getAnovaGrid(report);
  grok.shell.dockManager.dock(reportViewer.root, DG.DOCK_TYPE.FILL, hypoNode, 'Analysis');
}

/** Create dataframe with one-way ANOVA results. */
function getAnovaGrid(report: OneWayAnovaReport): DG.Grid {
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
    ['Source of variance', 'List of the explored variation sources'],
    ['SS', 'Sum of squares (SS)'],
    ['DF', 'Degrees of freedom (DF)'],
    ['MS', 'Mean square (MS)'],
    ['F', 'F-statistics (F)'],
    ['F-critical', `${report.significance}-critical value of F-statistics (F)`],
    ['p-value', `Probability to obtain F-statistics (F) greater than the actual observation.`],
  ]);

  grid.onCellTooltip(function(cell, x, y) {
    if (cell.isColHeader) {
      ui.tooltip.show(ui.divV([ui.p(tooltip.get(cell.tableColumn!.name)!)]), x, y);
      return true;
    }
  });

  grid.helpUrl = ANOVA_HELP_URL;

  return grid;
} // getOneWayAnovaDF

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
    onValueChanged: (col) => factor = col,
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
    onValueChanged: (col) => feature = col,
    filter: (col: DG.Column) => featureColNames.includes(col.name),
    nullable: false,
  });

  let significance = SIGNIFICANCE.DEFAULT;
  const signInput = ui.input.float('Alpha', {
    min: SIGNIFICANCE.MIN,
    max: SIGNIFICANCE.MAX,
    value: significance,
    nullable: false,
    tooltipText: 'Significance level',
    onValueChanged: (value) => {
      significance = value;
      runBtn.disabled = (significance <= SIGNIFICANCE.INFIMUM) || (significance >= SIGNIFICANCE.SUPREMUM);
    },
  });

  const dlg = ui.dialog({title: 'ANOVA', helpUrl: ANOVA_HELP_URL});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);
  dlg.addButton('Run', () => {
    dlg.close();

    try {
      const res = oneWayAnova(factor!, feature!, significance);
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

  dlg.add(factorInput)
    .add(featureInput)
    .add(signInput)
    .show();
} // runOneWayAnova
