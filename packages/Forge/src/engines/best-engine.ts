import * as DG from 'datagrok-api/dg';
import {Engine} from './engine';

const MANY_NUMERICAL_FEATURES = 5;

/** The built-in tool's suggestion: Chemprop for a molecule feature; XGBoost for a classification with categorical
 * and numerical features; PLS Regression for a regression with five or more numerical features, Linear Regression
 * for any other regression; XGBoost otherwise. A suggestion missing from [engines] gives the first engine. */
export function selectBestEngine(engines: Engine[], features: DG.Column[], target: DG.Column): Engine | undefined {
  const named = (name: string) => engines.find((e) => e.name === name) ?? engines[0];
  if (features.some((c) => c.semType === DG.SEMTYPE.MOLECULE))
    return named('Chemprop');
  const others = features.filter((c) => c.name !== target.name);
  const numerical = others.filter((c) => c.isNumerical).length;
  const categorical = others.filter((c) => c.isCategorical).length;
  const isRegression = target.isNumerical;
  if (!isRegression && categorical > 0 && numerical > 0)
    return named('XGBoost');
  if (isRegression)
    return named(numerical >= MANY_NUMERICAL_FEATURES ? 'PLS Regression' : 'Linear Regression');
  return named('XGBoost');
}
