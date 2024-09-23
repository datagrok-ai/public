import {expect} from '@datagrok-libraries/utils/src/test';
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const _A = 65;
const _Z = 90;


const isConjunction = (s: string) => s.toLowerCase().trim() === 'and' || s.toLowerCase().trim() === 'or';

const isEmpty = (s: string) => s === null || s === '';

const capitalize = (s: string) => isEmpty(s) ? s : (s[0].toUpperCase() + s.substring(1));

function isUpperCase(s: string, idx: number): boolean {
  const c = s.codePointAt(idx)!;
  return c >= _A && c <= _Z;
}

function splitCamelCase(s: string): string[] {
  const words: string[] = [];
  for (let start = 0, i = 1; i <= s.length; i++) {
    const lowerToUpper = i < s.length && i > 0 && isUpperCase(s, i) && !isUpperCase(s, i - 1);
    if (i === s.length || lowerToUpper || (isUpperCase(s, i) && i < s.length - 1 && !isUpperCase(s, i + 1))) {
      words[words.length] = s.substring(start, i);
      start = i;
    }
  }
  return words;
}

function camelCaseToWords(s: string, options?: {capitalizeFirst: boolean, capitalizeNext: boolean,
  capitalizeConjunctions: boolean, separator: string}): string | null {
  if (s == null)
    return null;
  if (s.toUpperCase() == s)
    return s;
  let sb = '';
  let prevWord;
  for (let word of splitCamelCase(s)) {
    word = (options?.capitalizeConjunctions ?? true) || !isConjunction(word) ? word : word.toLowerCase();
    if (sb.length == 0)
      sb += (options?.capitalizeFirst ?? true) ? capitalize(word) : word;
    else {
      if (!prevWord?.endsWith(options?.separator ?? ' '))
        sb += options?.separator ?? ' ';
      sb += (options?.capitalizeNext ?? false) ? capitalize(word) : word;
    }
    prevWord = word;
  }
  return sb;
}

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

export function customCaption(name: string, input: DG.InputBase, view: DG.View,
  customCaption: string, isInputForm: boolean = false): void {
  view.append(input.root);
  try {
    expect(input.caption, isInputForm ? customCaption : camelCaseToWords(customCaption));
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
      if (value === '')
        value = '[]';
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
