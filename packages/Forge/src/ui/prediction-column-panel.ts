import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {predictionTagOf} from '../apply/prediction-tag';
import {forgeDb} from '../generated/db';
import {CARD_COLUMNS, ForgeModelHandler, forgeModelHandler} from './model-handler';

export function isForgePrediction(col: DG.Column): boolean {
  return predictionTagOf(col) !== null;
}

/** **Predicted by**: the card of the model whose id the column's tag holds, or a note for a model of the earlier
 * predictive modeling tool. */
export function predictedByPane(col: DG.Column): HTMLElement {
  const tag = predictionTagOf(col);
  if (tag?.kind === 'legacy')
    return ui.divText(legacyPredictionText(tag.name));
  return ui.wait(async () => {
    const model = tag === null ? null :
      await forgeDb.models.query().where('id', '=', tag.id).select(...CARD_COLUMNS).first();
    return model === null ? ui.divText('The model that predicted this column is no longer available to you.') :
      forgeModelHandler.renderCard(ForgeModelHandler.rowOf(model));
  });
}

export function legacyPredictionText(name: string | null): string {
  return `Predicted by ${name === null ? '' : `"${name}", `}a model of the earlier predictive modeling tool. ` +
    'Forge cannot show it yet.';
}
