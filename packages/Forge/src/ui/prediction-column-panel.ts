import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PREDICTION_TAG} from '../constants';
import {forgeDb} from '../generated/db';
import {ForgeModelHandler, forgeModelHandler} from './model-handler';

// What the card, its tooltip and its commands (Download reads blob) read.
const CARD_COLUMNS = ['name', 'features', 'target_name', 'engine_name', 'task', 'row_count', 'blob'] as const;

export function isForgePrediction(col: DG.Column): boolean {
  return modelIdOf(col) !== '';
}

/** **Predicted by**: the card of the model whose id the column's tag holds. */
export function predictedByPane(col: DG.Column): HTMLElement {
  return ui.wait(async () => {
    const id = modelIdOf(col);
    const model = id === '' ? null :
      await forgeDb.models.query().where('id', '=', id).select(...CARD_COLUMNS).first();
    return model === null ? ui.divText('The model that predicted this column is no longer available to you.') :
      forgeModelHandler.renderCard(ForgeModelHandler.rowOf(model));
  });
}

function modelIdOf(col: DG.Column): string {
  return col.getTag(PREDICTION_TAG) ?? '';
}
