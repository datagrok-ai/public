import {expect} from '@datagrok-libraries/utils/src/test';
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function units(name: string, input: DG.InputBase, view: DG.View, selector: string, units: string): void {
  view.append(input.root);
  let value: string;
  try {
    value = (<HTMLInputElement>view.root.querySelector(selector)).innerText;
    expect(units, value);
  } catch (x) {
    throw new Error(name + ': ' + x);
  } finally {
    input.root.remove();
  }
}

export function caption(name: string, input: DG.InputBase, view: DG.View, selector: string): void {
  view.append(input.root);
  let value: string;
  try {
    value = (<HTMLInputElement>view.root.querySelector(selector)).innerText;
    expect(input.caption, value);
  } catch (x) {
    throw new Error(name + ': ' + x);
  } finally {
    input.root.remove();
  }
}

export function customCaption(name: string, input: DG.InputBase, view: DG.View, customCaption: string): void {
  view.append(input.root);
  try {
    expect(input.caption, customCaption);
  } catch (x) {
    throw new Error(name + ': ' + x);
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
      throw new Error(`"${name}": Element "${selector}" not found`);
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
    throw new Error(name + ': ' + x);
  } finally {
    input.enabled = true;
    input.root.remove();
  }
}

export function stringValue(name: string, input: DG.InputBase, selector: string, v: DG.View): void {
  v.root.innerHTML = '';
  v.append(input.root);
  let value: string;
  try {
    switch (name) {
    case 'multiChoiceInput':
      const node = v.root.querySelectorAll('.ui-input-multi-choice-checks>div');
      value = (<HTMLInputElement>v.root.querySelector('.ui-input-label')).innerText;
      for (let i = 0; i < node.length; i++) {
        const input = (<HTMLInputElement>node[i].querySelector('input'));
        if (input.checked)
          value = value.concat((<HTMLInputElement>node[i].querySelector('.ui-label')).innerText);
      }
      expect(input.stringValue, value);
      break;
    case 'columnInput':
      value = (<HTMLInputElement>v.root.querySelector('.d4-column-selector-column')).innerText;
      expect(input.stringValue, value);
      break;
    case 'columnsInput':
      value = (<HTMLInputElement>v.root.querySelector('.ui-input-column-names')).innerText.substring(3);
      expect(input.stringValue, value);
      break;
    case 'boolInput':
      value = (<HTMLInputElement>v.root.querySelector(`input${selector}`)).value;
      expect(input.stringValue, value);
      break;
    case 'switchInput':
      value = (<HTMLInputElement>v.root.querySelector(selector)?.querySelector('.ui-input-switch'))
        ?.classList.contains('ui-input-switch-on') ? 'true' : 'false';
      expect(input.stringValue, value);
      break;
    default:
      value = (<HTMLInputElement>v.root.querySelector(selector)).value;
      expect(input.stringValue, value);
      break;
    }
  } catch (x) {
    throw new Error(name + ': ' + x);
  } finally {
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
