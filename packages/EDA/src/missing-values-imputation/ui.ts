import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TITLE, KNN_IMPUTER, ERROR_MSG, HINT, MAX_INPUT_NAME_LENGTH} from './ui-constants';
import {SUPPORTED_COLUMN_TYPES, METRIC_TYPE, DISTANCE_TYPE, MetricInfo, DEFAULT, MIN_NEIGHBORS,
  impute, getMissingValsIndices, areThereFails, imputeFailed} from './knn-imputer';

/** Setting of the feature metric inputs */
type FeatureInputSettings = {
  defaultWeight: number,
  defaultMetric: METRIC_TYPE,
  availableMetrics: METRIC_TYPE[],
};

/** Return default setting of the feature metric inputs */
export function getFeatureInputSettings(type: DG.COLUMN_TYPE): FeatureInputSettings {
  switch (type) {
  case DG.COLUMN_TYPE.STRING:
  case DG.COLUMN_TYPE.DATE_TIME:
    return {
      defaultWeight: DEFAULT.WEIGHT,
      defaultMetric: METRIC_TYPE.ONE_HOT,
      availableMetrics: [METRIC_TYPE.ONE_HOT],
    };

  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT:
  case DG.COLUMN_TYPE.QNUM:
    return {
      defaultWeight: DEFAULT.WEIGHT,
      defaultMetric: METRIC_TYPE.DIFFERENCE,
      availableMetrics: [METRIC_TYPE.DIFFERENCE, METRIC_TYPE.ONE_HOT],
    };

  default:
    throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** Run the KNN missing values imputer */
export async function runKNNImputer(df?: DG.DataFrame): Promise<void> {
  /** current dataframe */
  df ??= grok.shell.t;

  if (df === null) {
    grok.shell.warning(ERROR_MSG.NO_DATAFRAME);
    return;
  }

  /** columns with missing values */
  const colsWithMissingVals = [] as DG.Column[];

  /** names of columns with missing values */
  const availableTargetColsNames = [] as string[];

  /** names of columns that can be used as features */
  const availableFeatureColsNames = [] as string[];

  // get columns with missing vals & available feature cols
  df.columns.toList()
    .filter((col) => SUPPORTED_COLUMN_TYPES.includes(col.type))
    .forEach((col) => {
      const misValsCount = col.stats.missingValueCount;
      if (misValsCount === col.length)
        return;

      availableFeatureColsNames.push(col.name);

      if (misValsCount > 0) {
        colsWithMissingVals.push(col);
        availableTargetColsNames.push(col.name);
      }
    });

  // get indices of missing values: col name -> array of indices
  const misValsInds = getMissingValsIndices(colsWithMissingVals);

  if (colsWithMissingVals.length === 0) {
    grok.shell.info(ERROR_MSG.NO_MISSING_VALUES);
    return;
  }

  if (availableFeatureColsNames.length === 1) {
    grok.shell.error(ERROR_MSG.ONE_AVAILABLE_FEATURE);
    return;
  }

  // In-place components
  let inPlace = DEFAULT.IN_PLACE > 0;
  const inPlaceInput = ui.input.bool(TITLE.IN_PLACE, {value: inPlace,
    onValueChanged: (value) => {inPlace = value ?? false;}});
  inPlaceInput.setTooltip(HINT.IN_PLACE);

  // Keep empty feature
  let keepEmpty = DEFAULT.KEEP_EMPTY > 0;
  const keepEmptyInput = ui.input.bool(TITLE.KEEP_EMPTY, {value: keepEmpty,
    onValueChanged: (value) => {keepEmpty = value ?? false;}});
  keepEmptyInput.setTooltip(HINT.KEEP_EMPTY);

  // Neighbors components
  let neighbors = DEFAULT.NEIGHBORS;
  const neighborsInput = ui.input.int(TITLE.NEIGHBORS, {
    value: neighbors,
    showPlusMinus: true,
    min: MIN_NEIGHBORS,
    nullable: false,
    onValueChanged: (value) => {
      if ((value !== null) && (value >= MIN_NEIGHBORS))
        neighbors = value;
      checkApplicability();
    },
  });
  neighborsInput.setTooltip(HINT.NEIGHBORS);

  // Distance components
  let distType = DISTANCE_TYPE.EUCLIDEAN;
  const distTypeInput: DG.ChoiceInput<DISTANCE_TYPE> = ui.input.choice(TITLE.DISTANCE, {
    value: distType,
    items: [DISTANCE_TYPE.EUCLIDEAN, DISTANCE_TYPE.MANHATTAN],
    onValueChanged: (value) => distType = value ?? DISTANCE_TYPE.EUCLIDEAN}) as DG.ChoiceInput<DISTANCE_TYPE>;
  distTypeInput.setTooltip(HINT.DISTANCE);

  // Target columns components (cols with missing values to be imputed)
  let targetColNames = colsWithMissingVals.map((col) => col.name);
  const targetColInput = ui.input.columns(TITLE.COLUMNS, {
    table: df,
    value: df.columns.byNames(availableTargetColsNames),
    onValueChanged: (value) => {
      targetColNames = value.map((col) => col.name);
      checkApplicability();
    },
    available: availableTargetColsNames,
  });
  targetColInput.setTooltip(HINT.TARGET);

  // Feature columns components
  let selectedFeatureColNames = availableFeatureColsNames as string[];
  const featuresInput = ui.input.columns(TITLE.FEATURES, {
    value: df.columns.byNames(availableFeatureColsNames),
    table: df, onValueChanged: (value) => {
      selectedFeatureColNames = value.map((col) => col.name);

      if (selectedFeatureColNames.length > 0) {
        checkApplicability();
        metricInfoInputs.forEach((div, name) => div.hidden = !selectedFeatureColNames.includes(name));
      } else
        hideWidgets();
    },
    available: availableFeatureColsNames,
  });
  featuresInput.setTooltip(HINT.FEATURES);

  /** Hide widgets (use if run is not applicable) */
  const hideWidgets = () => {
    dlg.getButton(TITLE.RUN).disabled = true;
    inPlaceInput.root.hidden = true;
    keepEmptyInput.root.hidden = true;
    neighborsInput.root.hidden = true;
    distDiv.hidden = true;
    metricsDiv.hidden = true;
  };

  /** Show widgets (use if run is applicable) */
  const showWidgets = () => {
    dlg.getButton(TITLE.RUN).disabled = (neighborsInput.value === null) || (neighborsInput.value < MIN_NEIGHBORS);
    distDiv.hidden = false;
    inPlaceInput.root.hidden = false;
    neighborsInput.root.hidden = false;
    distTypeInput.root.hidden = false;
    keepEmptyInput.root.hidden = !areThereFails(targetColNames, selectedFeatureColNames, misValsInds);
  };

  /** Check applicability of the imputation */
  const checkApplicability = () => {
    showWidgets();

    if (selectedFeatureColNames.length === 1) {
      targetColNames.forEach((name) => {
        if (selectedFeatureColNames[0] === name) {
          hideWidgets();
          grok.shell.warning(`${ERROR_MSG.ONE_FEATURE_SELECTED} the column '${name}'`);
        }
      });
    }

    if (targetColNames.length < 1)
      hideWidgets();
  };

  // Metrics components
  const featuresMetrics = new Map<string, MetricInfo>();
  const metricInfoInputs = new Map<string, HTMLElement>();
  const metricsDiv = ui.divV([]);
  metricsDiv.style.overflow = 'auto';

  // Create metrics UI
  availableFeatureColsNames.forEach((name) => {
    // initialization
    const type = df!.col(name)!.type as DG.COLUMN_TYPE;
    const settings = getFeatureInputSettings(type);
    featuresMetrics.set(name, {weight: settings.defaultWeight, type: settings.defaultMetric});

    // distance input
    const distTypeInput = ui.input.choice(name, {value: settings.defaultMetric,
      items: settings.availableMetrics, onValueChanged: (value) => {
        const distInfo = featuresMetrics.get(name) ?? {weight: settings.defaultWeight, type: settings.defaultMetric};
        distInfo.type = value ?? settings.defaultMetric;
        featuresMetrics.set(name, distInfo);
      }});
    distTypeInput.root.style.width = '50%';
    distTypeInput.setTooltip(HINT.METRIC);
    distTypeInput.root.hidden = true; // this input will be used further

    // The following should provide a slider (see th bug https://reddata.atlassian.net/browse/GROK-14431)
    const prop = DG.Property.fromOptions({
      'name': name.length < MAX_INPUT_NAME_LENGTH ? name : name.slice(0, MAX_INPUT_NAME_LENGTH).concat('...'),
      'inputType': 'Float',
      'min': 0,
      'max': 10,
      // @ts-ignore
      'showSlider': true,
      'step': 1,
    });
    const weightInput = ui.input.forProperty(prop);
    weightInput.value = settings.defaultWeight;
    weightInput.onChanged.subscribe((value) => {
      const distInfo = featuresMetrics.get(name) ?? {weight: settings.defaultWeight, type: settings.defaultMetric};
      distInfo.weight = value ?? settings.defaultWeight;
      featuresMetrics.set(name, distInfo);
    });
    ui.tooltip.bind(weightInput.captionLabel, name);
    ui.tooltip.bind(weightInput.input, HINT.WEIGHT);

    metricInfoInputs.set(name, weightInput.root);
    metricsDiv.append(weightInput.root);
  });

  // The main dialog
  const dlg = ui.dialog({title: TITLE.KNN_IMPUTER, helpUrl: KNN_IMPUTER});
  grok.shell.v.root.appendChild(dlg.root);

  metricsDiv.hidden = true;
  keepEmptyInput.root.hidden = !areThereFails(targetColNames, selectedFeatureColNames, misValsInds);

  // Icon showing/hiding metrics UI
  const settingsIcon = ui.icons.settings(() => {metricsDiv.hidden = !metricsDiv.hidden;}, HINT.METRIC_SETTINGS);

  const distDiv = ui.divH([distTypeInput.root, settingsIcon]);

  let resolve: (value: void | PromiseLike<void>) => void;
  let reject: (reason?: any) => void;
  let okClicked = false;
  const promise = new Promise<void>((res, rej) => {
    resolve = res;
    reject = rej;
  });

  dlg.addButton(TITLE.RUN, () => {
    okClicked = true;
    dlg.close();
    availableFeatureColsNames.filter((name) => !selectedFeatureColNames.includes(name))
      .forEach((name) => featuresMetrics.delete(name));

    try {
      const failedToImpute = impute(df!, targetColNames, featuresMetrics, misValsInds, distType, neighbors, inPlace);

      if (!keepEmpty)
        imputeFailed(df!, failedToImpute);
      resolve();
    } catch (err) {
      if (err instanceof Error)
        grok.shell.error(`${ERROR_MSG.KNN_FAILS}: ${err.message}`);
      else
        grok.shell.error(`${ERROR_MSG.KNN_FAILS}: ${ERROR_MSG.CORE_ISSUE}`);
      reject(err);
    }
  }, undefined, HINT.RUN);

  dlg.add(targetColInput)
    .add(featuresInput)
    .add(distDiv)
    .add(metricsDiv)
    .add(neighborsInput)
    .add(inPlaceInput)
    .add(keepEmptyInput)
    .show()
    .onClose.subscribe(() => !okClicked && resolve());

  return promise;
} // runKNNImputer
