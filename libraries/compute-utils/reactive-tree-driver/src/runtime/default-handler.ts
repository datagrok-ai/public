import {RestrictionType} from '../data/common-types';
import {LinkController} from './LinkControllers';

export function defaultLinkHandler(
  ctrlInstance: LinkController,
  inputs: string[],
  outputs: string[],
  defaultRestrictions?: Record<string, RestrictionType>,
) {
  for (let i = 0; i < Math.min(inputs.length, outputs.length); i++) {
    const input = ctrlInstance.getFirst(inputs[i]);
    const restriction = defaultRestrictions?.[outputs[i]];
    ctrlInstance.setAll(outputs[i], input, restriction);
  }
  ctrlInstance.close();
}
