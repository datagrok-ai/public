// Analysis of Variances (ANOVA): UI

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {oneWayAnova, OneWayAnovaReport} from './anova-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];

const ANOVA_HELP_URL = '/help/explore/anova';

/** Add one-way ANOVA results */
function addVizualization(df: DG.DataFrame, factors: DG.Column, values: DG.Column, report: OneWayAnovaReport): void {
  const test = report.anovaTable.fStat > report.fCritical;

  const shortConclusion = test ?
    `"${factors.name}" affects the "${values.name}"` :
    `"${factors.name}" doesn't affect the "${values.name}"`;

  const view = grok.shell.getTableView(df.name);
  const boxPlot = DG.Viewer.boxPlot(df, {
    categoryColumnName: factors.name,
    valueColumnName: values.name,
    showPValue: false,
    showStatistics: false,
    description: shortConclusion,
  });
  const boxPlotNode = view.dockManager.dock(boxPlot.root, DG.DOCK_TYPE.RIGHT, null, 'ANOVA');

  const hypoMd = ui.markdown(`**H0:** the "${factors.name}" 
  factor does not produce a significant difference in the "${values.name}" feature.`);
  ui.tooltip.bind(hypoMd, 'Null hypothesis');

  const testMd = ui.markdown(`**Test result:** ${test ?
    'there is significant differences among sample averages.' :
    'there is no significant differences among sample averages.'}`,
  );

  const sign = test ? '>' : '<';
  ui.tooltip.bind(testMd, () => ui.markdown(`${test ?
    'Reject the null hypothesis' :
    'Fail to reject the null hypothesis'}, since F ${sign} F-critical: 
    ${report.anovaTable.fStat.toFixed(2)} ${sign} ${report.fCritical.toFixed(2)}.`));

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

  const factorColNames = [] as string[];
  const featureColNames = [] as string[];

  for (const col of df.columns) {
    if (col.stats.missingValueCount < 1) {
      if (FEATURE_TYPES.includes(col.type))
        featureColNames.push(col.name);
      else if (FACTOR_TYPES.includes(col.type))
        factorColNames.push(col.name);
    }
  }

  if (factorColNames.length < 1) {
    grok.shell.warning(ui.markdown(`No acceptable factor columns:

    - no missing values
    - type: ${FACTOR_TYPES.join(', ')} 
    - at least two categories`,
    ));
    return;
  }

  if (featureColNames.length < 1) {
    grok.shell.warning(ui.markdown(`No acceptable feature columns:
    
    - no missing values
    - type: ${FEATURE_TYPES.join(', ')}`,
    ));
    return;
  }

  let factor = df.col(factorColNames[0]);
  const factorInput = ui.columnInput('Category', df, factor, () => factor = factorInput.value, {
    filter: (col: DG.Column) => factorColNames.includes(col.name),
  });
  factorInput.setTooltip('Column with factor values');

  let feature = df.col(featureColNames[0]);
  const featureInput = ui.columnInput('Feature', df, feature, () => feature = featureInput.value, {
    filter: (col: DG.Column) => featureColNames.includes(col.name),
  });
  featureInput.setTooltip('Column with feature values');

  let significance = 0.05;
  const signInput = ui.input.forProperty(DG.Property.fromOptions({
    name: 'significance',
    defaultValue: significance,
    inputType: 'Float',
    min: 0.01,
    max: 0.99,
  }));
  signInput.onChanged(() => {
    significance = signInput.value;
    runBtn.disabled = (significance <= 0) || (significance >= 1);
  });
  signInput.setTooltip('Specifies the criterion used for rejecting the null hypothesis.');

  let validate = true;
  const validateInput = ui.input.forProperty(DG.Property.fromOptions({
    name: 'validate',
    defaultValue: validate,
    inputType: 'Bool',
  }));
  validateInput.onChanged(() => validate = validateInput.value);
  validateInput.setTooltip('Indicates whether to check applicability of ANOVA.');

  const dlg = ui.dialog({title: 'ANOVA', helpUrl: ANOVA_HELP_URL});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);
  dlg.addButton('Run', () => {
    dlg.close();

    try {
      const res = oneWayAnova(factor!, feature!, significance, validate);
      addVizualization(df, factor!, feature!, res);
    } catch (error) {
      if (error instanceof Error) {
        grok.shell.warning(getWarning(error.message));

        view.addViewer(DG.VIEWER.BOX_PLOT, {
          categoryColumnName: factor!.name,
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
    .add(validateInput)
    .show();
} // runOneWayAnova
