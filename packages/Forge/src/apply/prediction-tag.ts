import * as DG from 'datagrok-api/dg';
import {PREDICTION_TAG} from '../constants';

export type PredictionTag = {kind: 'forge', id: string} | {kind: 'legacy', name: string | null};

const LEGACY_MARKUP = /#\{x\.[0-9a-f-]{36}\."(.*)"\}/;

/** The model a prediction column names in {@link PREDICTION_TAG}: a Forge model's id, or the markup an earlier
 * predictive modeling tool's model wrote there; null when the column has no such tag. */
export function predictionTagOf(col: DG.Column): PredictionTag | null {
  const value = col.getTag(PREDICTION_TAG) ?? '';
  if (value === '')
    return null;
  const legacy = LEGACY_MARKUP.exec(value);
  return legacy === null ? {kind: 'forge', id: value} : {kind: 'legacy', name: legacy[1] === '' ? null : legacy[1]};
}
