import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


export function caption(name: string, input: DG.InputBase, view: DG.View, selector: string): void {
  view.append(input.root);
  let value: string;
  try {
    value = (<HTMLInputElement>view.root.querySelector(selector)).innerText;
    expect(input.caption, value);
  } catch (x) {
    throw name + ': ' + x;
  } finally {
    input.root.remove();
  }
}

export function checkHTMLElement(name: string, root: HTMLElement, v: DG.View, selectors: string | string[]): void {
  if (typeof selectors === 'string')
    selectors = [selectors];

  v.append(root);
  selectors.forEach((selector) => {
    const e = v.root.querySelector(selector);
    if (e == undefined)
      throw `"${name}": Element "${selector}" not found`;
  });
  root.remove();
}

export function enabled(name: string, input: DG.InputBase, v: DG.View, selector: string): void {
  v.append(input.root);
  let value: boolean = true;
  input.enabled = false;
  try {
    if (<HTMLInputElement>v.root.querySelector(selector))
      value = false;
    expect(input.enabled, value);
  } catch (x) {
    throw name + ': ' + x;
  } finally {
    input.enabled = true;
    input.root.remove();
  }
}
