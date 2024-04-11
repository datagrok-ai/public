// Analysis of Variances (ANOVA): UI

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {oneWayAnova} from './anova-tools';

const FEATURE_TYPES = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as string[];
const FACTOR_TYPES = [DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.BOOL] as string[];


/** Add one-way ANOVA results */
function addOneWayAnovaVizualization(
  table: DG.DataFrame, factors: DG.Column, values: DG.Column, anova: DG.DataFrame,
): void {
  const view = grok.shell.getTableView(table.name);
  view.addViewer(DG.VIEWER.BOX_PLOT, {
    categoryColumnName: factors.name,
    valueColumnName: values.name,
    showPValue: false,
    showValueSelector: false,
    showCategorySelector: false,
    showStatistics: false,
  });
  view.addViewer(DG.Viewer.grid(anova));
}

/** Return warning div */
function getWarning(msg: string): HTMLElement {
  return ui.divV([
    ui.markdown(`ANOVA cannot be performed:
        
    ${msg}`),
    ui.link('Learn more',
      () => window.open('https://en.wikipedia.org/wiki/Analysis_of_variance#Assumptions', '_blank'),
      'Open link in a new tab',
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

  const dlg = ui.dialog({title: 'ANOVA', helpUrl: '/help/explore/anova'});
  const view = grok.shell.getTableView(df.name);
  view.root.appendChild(dlg.root);
  dlg.addButton('Run', () => {
    dlg.close();

    try {
      const res = oneWayAnova(factor!, feature!, significance, validate);
      addOneWayAnovaVizualization(df, factor!, feature!, res);
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
