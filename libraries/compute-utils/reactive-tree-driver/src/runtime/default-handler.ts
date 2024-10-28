import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
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
    if (input instanceof DG.DataFrame)
      ctrlInstance.setAll(outputs[i], input.clone(), restriction);
    else
      ctrlInstance.setAll(outputs[i], input, restriction);
  }
  ctrlInstance.close();
}
