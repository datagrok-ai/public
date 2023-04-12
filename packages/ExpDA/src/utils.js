// utils.js

// Inputs correctness check tools

//Limitation constants
const COMP_MIN = 1;
const SAMPLES_COUNT_MIN = 1;
const FEATURES_COUNT_MIN = 1;
const PERCENTAGE_MIN = 0;
const PERCENTAGE_MAX = 100;

// Error messages
const COMP_POSITVE_MES = 'components must be positive.';
const COMP_EXCESS = 'components must not be greater than feautures count.';
const INCORERRECT_MIN_MAX_MES = 'min must be less than max.';
const INCORERRECT_FEATURES_MES = 'features must be positive.';
const INCORERRECT_SAMPLES_MES = 'samples must be positive.';
const INCORERRECT_PERCENTAGE_MES = 'violators percentage must be from the range from 0 to 100.';

// Check components count (PCA, PLS)
export function checkComponenets(features, components) {
  if (components < COMP_MIN)
    throw new Error(COMP_POSITVE_MES);

  if (components > features.length)
    throw new Error(COMP_EXCESS);
}

// Check whether data processing requires much time (PCA, PLS)
export function isProcExpensive(table, features, max) {
  return ((table.rowCount > max) || (features.length > max));
}

// Check inputs of data for SVM testing generator
export function checkGeneratorSVMinputs(samplesCount, featuresCount, min, max, violatorsPercentage) 
{
  if (min >= max)
    throw new Error(INCORERRECT_MIN_MAX_MES);
  
  if (featuresCount < FEATURES_COUNT_MIN) 
    throw new Error(INCORERRECT_FEATURES_MES);

  if (samplesCount < SAMPLES_COUNT_MIN) 
    throw new Error(INCORERRECT_SAMPLES_MES);

  if ((violatorsPercentage < PERCENTAGE_MIN) || (violatorsPercentage > PERCENTAGE_MAX))
    throw new Error(INCORERRECT_PERCENTAGE_MES);
}

// Check whether SVM computation requires much time
export function isDataBigForSVM(featureCols, featuresMax, labelsCol, labelsMax) {
  if (featureCols.length > featuresMax)
    return true;

  if (labelsCol.length > labelsMax)
    return true;
}
