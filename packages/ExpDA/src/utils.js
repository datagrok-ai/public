// inputChecks.js

// Dataframe size limitation constants
const COMP_MIN = 1;

// Error messages
const COMP_POSITVE_MES = 'components must be positive.';
const COMP_EXCESS = 'components must not be greater than feautures count.';

// Check components count (PCA, PLS)
export function checkComponenets(features, components) {
  if (components < COMP_MIN)
    throw new Error(COMP_POSITVE_MES);

  if (components > features.length)
    throw new Error(COMP_EXCESS);
}

// Check whether data processing requires much time
export function isProcExpensive(table, features, max) {
  return ((table.rowCount > max) || (features.length > max));
}
