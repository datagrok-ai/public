import * as DG from 'datagrok-api/dg';
import {waitUntil} from '@datagrok-libraries/test/src/test';
import type {
  ConsistencyInfo, RestrictionType, ValidationIconInput,
} from '@datagrok-libraries/webcomponents/src/ValidationIcon/ValidationIcon';

export const scriptWith = (...inputs: string[]) => DG.Script.create([
  '//name: webcomponentsFixture', '//language: javascript', ...inputs.map((input) => `//input: ${input}`),
  '//output: double out',
].join('\n'));

export function typeInto(input: DG.InputBase, value: string) {
  const el = input.input as HTMLInputElement;
  el.focus();
  Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!.call(el, value);
  el.dispatchEvent(new Event('input', {bubbles: true}));
  el.dispatchEvent(new Event('change', {bubbles: true}));
  el.blur();
}

export const consistency = (restriction: RestrictionType, inconsistent = false, assignedValue: any = 1):
  ConsistencyInfo => ({restriction, inconsistent, assignedValue});

export const ADVICE: ValidationIconInput = {validation: {warnings: [{
  description: 'Too high',
  actions: [{actionName: 'Set to 10', action: 'action-1', additionalParams: {val: 10}}],
}]}};

export const iconOf = (root: HTMLElement) => root.querySelector('dg-validation-icon i') as HTMLElement | null;

export const findLink = (root: ParentNode, text: string) => Array.from(root.querySelectorAll('a.ui-link'))
  .find((el) => el.textContent?.trim() === text) as HTMLElement | undefined;

export async function openPopover(root: HTMLElement) {
  await waitUntil(() => iconOf(root) != null, 'validation icon not rendered', 3000);
  iconOf(root)!.click();
  await waitUntil(() => !!root.querySelector('.dg-validation-icon-popover')?.childElementCount,
    'validation popover not rendered', 3000);
  return root.querySelector('.dg-validation-icon-popover') as HTMLElement;
}
