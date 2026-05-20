import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RestrictionType} from '../data/common-types';
import {LinkController} from './LinkControllers';
import {TemplateId} from '../RuntimeControllers';

export type Slot = {kind: 'bare', name: string} | {kind: 'template', name: TemplateId};

export function defaultLinkHandler(
  ctrlInstance: LinkController,
  inputs: Slot[],
  outputs: Slot[],
  defaultRestrictions?: Record<string, RestrictionType> | RestrictionType,
) {
  for (let i = 0; i < Math.min(inputs.length, outputs.length); i++) {
    const inS = inputs[i];
    const outS = outputs[i];
    if (inS.kind === 'template' && outS.kind === 'template') {
      ctrlInstance.propagateTemplatePair(inS.name, outS.name, defaultRestrictions);
    } else if (inS.kind === 'bare' && outS.kind === 'bare') {
      if (ctrlInstance.callInputs.has(inS.name))
        continue;
      const input = ctrlInstance.getFirst(inS.name);
      const restriction = typeof defaultRestrictions === 'string' ?
        defaultRestrictions :
        (defaultRestrictions?.[outS.name] ?? defaultRestrictions?.['*']);
      ctrlInstance.setAll(outS.name, input, restriction);
    } else {
      throw new Error(`Link ${ctrlInstance.id}: slot ${i} mismatches — ${inS.kind} on input vs ${outS.kind} on output. Template-to-bare propagation is not allowed.`);
    }
  }
  ctrlInstance.close();
}
