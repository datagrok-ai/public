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

export function waitForHTMLCollection(selector: string, wait=3000): Promise<HTMLCollection> {
  return new Promise((resolve, reject) => {
    if (document.querySelector(selector) !== null) {
      if (document.querySelector(selector)!.children.length !== 0)
        return resolve(document.querySelector(selector)!.children);
    }

    const observer = new MutationObserver(() => {
      if (document.querySelector(selector) !== null) {
        if (document.querySelector(selector)!.children.length !== 0) {
          clearTimeout(timeout);
          observer.disconnect();
          resolve(document.querySelector(selector)!.children);
        }
      }
    });

    const timeout = setTimeout(() => {
      observer.disconnect();
      reject(new Error(`cannot find ${selector}!`));
    }, wait,
    );

    observer.observe(document.body, {
      childList: true,
      subtree: true,
    });
  });
}

export function waitForHTMLElement(selector: string, regex: RegExp, error: string, wait=3000): Promise<HTMLElement> {
  return new Promise((resolve, reject) => {
    if (document.querySelector(selector) !== null) {
      if (regex.test((document.querySelector(selector) as HTMLElement).innerText))
        return resolve(document.querySelector(selector) as HTMLElement);
    }

    const observer = new MutationObserver(() => {
      if (document.querySelector(selector) !== null) {
        if (regex.test((document.querySelector(selector) as HTMLElement).innerText)) {
          clearTimeout(timeout);
          observer.disconnect();
          resolve(document.querySelector(selector) as HTMLElement);
        }
      }
    });

    const timeout = setTimeout(() => {
      observer.disconnect();
      reject(new Error(error));
    }, wait,
    );

    observer.observe(document.body, {
      childList: true,
      subtree: true,
    });
  });
}
