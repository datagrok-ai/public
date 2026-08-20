import {span} from '../core/elements.js';
import {DynamicInput, DynamicInputOptions} from '../components/inputs/dynamic-input.js';
import {ChipOptions, EntityCard, entityCard} from './entity.js';

export interface MetaInputOptions extends Omit<DynamicInputOptions<any>, 'render'>, ChipOptions {}

/** A read-only input showing whatever object it holds as the platform's handler card — the
 * counterpart of Dart's `MetaInput` (`dynamic_input.dart:44-53`, `EntityMeta.renderPreview`), down
 * to the `<null>` marker for an empty value. Lives in `u2/dg` because the card comes from the
 * ObjectHandler registry.
 *
 * The card is a component, so each swap disposes the previous one: nothing owns what a renderer
 * builds inside an effect, since `Scope.ambient` is only set while a tree is being built. */
export function metaInput(options: MetaInputOptions = {}): DynamicInput<any> {
  let card: EntityCard | null = null;
  const input = new DynamicInput<any>({
    ...options,
    render: (value) => {
      card?.dispose();
      card = value == null ? null : entityCard(value, options);
      return card === null ? span('<null>', 'u2-input-null') : card.root;
    },
  });
  input.own(() => card?.dispose());
  return input;
}
