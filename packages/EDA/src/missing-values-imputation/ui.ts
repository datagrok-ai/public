import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { TITLE, LINK, ERROR_MSG, HINT } from './ui-constants';
import { DISTANCE_TYPE, DistanceInfo, DEFAULT, MIN_NEIGHBORS, impute } from "./knn-imputer";

/** */
const SUPPORTED_COLUMN_TYPES = [
  DG.COLUMN_TYPE.INT,
  DG.COLUMN_TYPE.FLOAT,
  DG.COLUMN_TYPE.STRING,
  DG.COLUMN_TYPE.DATE_TIME,
  DG.COLUMN_TYPE.QNUM,
];

/** */
type FeatureInputSettings = {
  defaultWeight: number,
  defaultDistance: DISTANCE_TYPE,
  availableDistances: DISTANCE_TYPE[], 
};

/** */
function getFeatureInputSettings(type: DG.COLUMN_TYPE): FeatureInputSettings {
  switch (type) {
    case DG.COLUMN_TYPE.STRING:
    case DG.COLUMN_TYPE.DATE_TIME:
      return {
        defaultWeight: DEFAULT.WEIGHT,
        defaultDistance: DISTANCE_TYPE.ONE_HOT,
        availableDistances: [DISTANCE_TYPE.ONE_HOT]
      };

    case DG.COLUMN_TYPE.INT:
    case DG.COLUMN_TYPE.FLOAT:
    case DG.COLUMN_TYPE.QNUM:
      return {
        defaultWeight: DEFAULT.WEIGHT,
        defaultDistance: DISTANCE_TYPE.EUCLIDEAN,
        availableDistances: [DISTANCE_TYPE.EUCLIDEAN, DISTANCE_TYPE.MANHATTAN, DISTANCE_TYPE.ONE_HOT]
      };

    default:
      throw new Error(ERROR_MSG.UNSUPPORTED_COLUMN_TYPE);
  }
}

/** */
function isSupportedColumnType(col: DG.Column): boolean {
  return SUPPORTED_COLUMN_TYPES.includes(col.type as DG.COLUMN_TYPE);
}

/** Check whether missing values imputation can be applied */
function canImputationBeApplied(col: DG.Column): boolean {
  return isSupportedColumnType(col) && (col.stats.missingValueCount > 0);
}

/** */
export function runKNNImputer(): void {
  let df: DG.DataFrame | null = grok.shell.t;

  if (df === null)
    throw new Error(ERROR_MSG.NO_DATAFRAME);

  const colsWithMissingVals = [] as DG.Column[];
  const availableColsNames = [] as string[];

  df.columns.toList()
    .filter((col) => isSupportedColumnType(col))
    .forEach((col) => {
      availableColsNames.push(col.name);

      if (canImputationBeApplied(col))
        colsWithMissingVals.push(col);
    });

  if (colsWithMissingVals.length === 0) {
    grok.shell.info(ERROR_MSG.NO_MISSING_VALUES);
    return;
  }

  if (colsWithMissingVals.length === 1) {
    grok.shell.error(ERROR_MSG.KNN_CANNOT_BE_APPLIED);
    return;
  }

  let inPlace = DEFAULT.IN_PLACE > 0;
  const inPlaceInput = ui.boolInput(TITLE.IN_PLACE, inPlace, () => { inPlace = inPlaceInput.value ?? false;});
  inPlaceInput.setTooltip(HINT.IN_PLACE);

  let neighbors = DEFAULT.NEIGHBORS;
  const neighborsInput = ui.intInput(TITLE.NEIGHBORS, neighbors, () => {
    const val = neighborsInput.value;
    if (val === null)
      neighborsInput.value = neighbors;
    else if (val > MIN_NEIGHBORS)
      neighbors = val;
  });
  neighborsInput.setTooltip(HINT.NEIGHBORS);

  let targetCol: DG.Column | null = null;
  const targetColInput = ui.columnInput(TITLE.COLUMN, df, null, () => {
    hideWidgets();
    targetCol = targetColInput.value;

    if (targetCol !== null)
      featuresInput.root.hidden = false;
  }, {filter: (col: DG.Column) => colsWithMissingVals.includes(col)});
  targetColInput.setTooltip(HINT.FEATURES);

  let selectedFeatureColNames = [] as string[];
  const featuresInput = ui.columnsInput(TITLE.FEATURES, df, () => {
    hideWidgets();
    featuresInput.root.hidden = false;
    selectedFeatureColNames = featuresInput.value.map((col) => col.name).filter((name) => name !== targetCol!.name);
    featuresInput.value = featuresInput.value.filter((col) => col.name !== targetCol!.name);

    if (selectedFeatureColNames.length > 0) {
      dlg.getButton(TITLE.RUN).hidden = false;
      inPlaceInput.root.hidden = false;
      neighborsInput.root.hidden = false;

      selectedFeatureColNames.forEach((name) => {
        const uiElem = distInfoInputs.get(name);
        if (uiElem !== undefined)
          uiElem.hidden = false;
      });
    }
  }, {available: availableColsNames});
  featuresInput.setTooltip(HINT.FEATURES);

  const featuresDist = new Map<string, DistanceInfo>();
  const distInfoInputs = new Map<string, HTMLDivElement>();
  const metricSpecificationInputs = ui.divV([]);

  availableColsNames.forEach((name) => {    
    const type = df!.col(name)!.type as DG.COLUMN_TYPE;
    const settings = getFeatureInputSettings(type);
    featuresDist.set(name, {weight: settings.defaultWeight, type: settings.defaultDistance});

    const distTypeInput = ui.choiceInput(name, settings.defaultDistance, settings.availableDistances, () => {
      const distInfo = featuresDist.get(name) ?? {weight: settings.defaultWeight, type: settings.defaultDistance};
      distInfo.type = distTypeInput.value ?? settings.defaultDistance;
      featuresDist.set(name, distInfo);
    });
    distTypeInput.root.style.width = '50%';
    distTypeInput.setTooltip(HINT.DISTANCE);

    const weightInput = ui.floatInput(TITLE.WEIGHT, settings.defaultWeight, () => {
      const distInfo = featuresDist.get(name) ?? {weight: settings.defaultWeight, type: settings.defaultDistance};
      distInfo.weight = weightInput.value ?? settings.defaultWeight;
      featuresDist.set(name, distInfo);
    });
    weightInput.root.style.width = '10%';
    weightInput.setTooltip(HINT.WEIGHT);
    
    const div = ui.divH([distTypeInput.root, weightInput.root]);
    distInfoInputs.set(name, div);
    metricSpecificationInputs.append(div);
  });  

  const dlg = ui.dialog({title: TITLE.KNN_IMPUTER, helpUrl: LINK.KNN_IMPUTER});
  grok.shell.v.root.appendChild(dlg.root);

  dlg.addButton(TITLE.RUN, () => {
      dlg.close();

      availableColsNames
        .filter((name) => !selectedFeatureColNames.includes(name))
        .forEach((name) => featuresDist.delete(name));

      impute(df!, targetCol!, featuresDist, neighbors, inPlace);
    })
    .add(targetColInput)
    .add(featuresInput)
    .add(metricSpecificationInputs)
    .add(neighborsInput)
    .add(inPlaceInput)
    .show();
  
  const hideWidgets = () => {
    dlg.getButton(TITLE.RUN).hidden = true;
    inPlaceInput.root.hidden = true;
    neighborsInput.root.hidden = true;
    featuresInput.root.hidden = true;
    distInfoInputs.forEach((div) => div.hidden = true);
  };

  hideWidgets();
}
