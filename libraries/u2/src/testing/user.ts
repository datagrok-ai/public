import {Control} from 'datagrok-api/u2core';
import {Names, Query} from './query.js';
import {Reach} from './probe.js';

/** What a gesture acts on: an element or a control (its root). */
export type Target = Element | {root: Element};

/** Options of {@link UserGestures.type}. */
export interface TypeOptions {
  /** Selects what the field holds first, so the text replaces it. */
  replace?: boolean;
}

/** The gestures of {@link user}. Each first checks that a user could make it — the element is on the page, shown,
 * not inert or disabled, takes the pointer, and is not covered where the page has layout (`DG.U2.operability`) — and
 * throws what stops them; then it dispatches what the browser would. */
export interface UserGestures {
  /** Presses the pointer on the element and lets it go: pointer and mouse events, focus, `click`. */
  click(target: Target, init?: MouseEventInit): void;
  /** Clicks into a text field (the element, or the field of the input it is) and types `text` key by key: `keydown`,
   * the text at the caret, `input`, `keyup`. A read-only or disabled field refuses it. */
  type(target: Target, text: string, options?: TypeOptions): void;
  /** Presses `key` on the element, focused first: `keydown`, `keyup`; Enter or Space on a button clicks it. */
  press(target: Target, key: string, init?: KeyboardEventInit): void;
  /** Picks the item labelled `item`: an option of the element's `<select>`, else an option, tree item, menu item or
   * tab inside it or in a popup it opened, clicked. Throws listing what it offers when nothing is labelled so. */
  pick(target: Target, item: string): void;
}

const FOCUSABLE = 'a[href], button, input, select, textarea, summary, [tabindex], [contenteditable="true"]';
const TEXT_FIELD = 'input:not([type="checkbox"]):not([type="radio"]), textarea, [contenteditable="true"]';
const ITEMS = '[role="option"], [role="treeitem"], [role="menuitem"], [role="menuitemcheckbox"], ' +
  '[role="menuitemradio"], [role="tab"]';

/** A user's gestures, refused where a user would be refused: `user.click(byName(dialog, 'Edit Campaign Name'))`
 * throws `button "Edit Campaign Name" takes no pointer: pointer-events: none on div.u2-input-options` where no user
 * could press it, although a synthetic event would still reach its handler. The same code runs in Node and in a
 * browser package test.
 * @example
 * user.click(byName(dialog, 'Edit Campaign Name'));
 * user.type(byName(dialog, 'Name'), 'KHL-4', {replace: true});
 * user.press(byName(dialog, 'Name'), 'Enter');
 * user.pick(byName(dialog, 'submit-function'), 'Demo File Submit'); */
export const user: UserGestures = {
  click(target: Target, init?: MouseEventInit): void {
    Gesture.click(Query.element(target), 'click', init);
  },

  type(target: Target, text: string, options?: TypeOptions): void {
    const el = Query.element(target);
    const field = el.matches(TEXT_FIELD) ? el : el.querySelector<HTMLElement>(TEXT_FIELD);
    if (field === null)
      throw new Error(`user.type: ${Names.called(el)} has no text field to type into`);
    const focused = document.activeElement === field;
    if (focused)
      Gesture.reach(field, 'type', false);
    else
      Gesture.click(field, 'type');
    const input = field as HTMLInputElement;
    if (input.readOnly || input.disabled)
      throw new Error(`user.type: ${Names.called(field)} is read-only: a user cannot type into it`);
    const caret = Caret.of(input, options?.replace === true, focused);
    for (const ch of text) {
      const key = {key: ch, bubbles: true, cancelable: true, composed: true};
      if (field.dispatchEvent(new KeyboardEvent('keydown', key))) {
        caret.insert(ch);
        field.dispatchEvent(typeof InputEvent === 'function' ?
          new InputEvent('input', {bubbles: true, composed: true, data: ch, inputType: 'insertText'}) :
          new Event('input', {bubbles: true, composed: true}));
      }
      field.dispatchEvent(new KeyboardEvent('keyup', key));
    }
  },

  press(target: Target, key: string, init?: KeyboardEventInit): void {
    const el = Gesture.reach(Query.element(target), 'press', false);
    if (document.activeElement !== el && el.matches(FOCUSABLE))
      el.focus();
    const opts = {key, bubbles: true, cancelable: true, composed: true, ...init};
    const go = el.dispatchEvent(new KeyboardEvent('keydown', opts));
    const activates = key === 'Enter' ? el.matches('button, a[href]') :
      key === ' ' && el.matches('button, input[type="checkbox"], input[type="radio"]');
    if (go && activates)
      el.click();
    el.dispatchEvent(new KeyboardEvent('keyup', opts));
  },

  pick(target: Target, item: string): void {
    const el = Query.element(target);
    const select = el.matches('select') ? el as HTMLSelectElement : el.querySelector('select');
    if (select !== null) {
      const options = Array.from(select.querySelectorAll('option'));
      const option = options.find((o) => o.value === item || Names.same(o.textContent ?? '', item));
      if (option === undefined) {
        throw new Error(`user.pick: ${Names.called(select)} offers no "${item}"; it offers ` +
          `${options.map((o) => `"${Names.squash(o.textContent ?? '')}"`).join(', ') || 'nothing'}`);
      }
      if (option.disabled)
        throw new Error(`user.pick: "${item}" is disabled in ${Names.called(select)}`);
      Gesture.click(select, 'pick');
      select.value = option.value;
      option.selected = true;
      select.dispatchEvent(new Event('input', {bubbles: true, composed: true}));
      select.dispatchEvent(new Event('change', {bubbles: true}));
      return;
    }
    const owner = Control.forElement(el)?.root.id;
    const popups = owner ? Array.from(document.querySelectorAll('[data-u2-owner-id]'))
      .filter((p) => p.getAttribute('data-u2-owner-id') === owner) : [];
    const items = [el, ...popups].flatMap((root) => Query.all(root, ITEMS));
    const found = items.filter((i) => Names.same(Names.accessible(i), item));
    if (found.length !== 1) {
      throw new Error(`user.pick: ${found.length === 0 ? 'nothing' : `${found.length} items`} in ` +
        `${Names.called(el)} ${found.length === 0 ? 'is' : 'are'} labelled "${item}"; it offers ` +
        `${[...new Set(items.map((i) => `"${Names.accessible(i)}"`))].join(', ') || 'nothing'}`);
    }
    Gesture.click(found[0], 'pick');
  },
};

class Gesture {
  static reach(el: HTMLElement, gesture: string, pointer: boolean): HTMLElement {
    const why = Reach.blocked(el, pointer);
    if (why !== undefined)
      throw new Error(`user.${gesture}: ${why}`);
    return el;
  }

  static click(target: HTMLElement, gesture: string, init?: MouseEventInit): void {
    const el = Gesture.reach(target, gesture, true);
    const r = el.getBoundingClientRect();
    const base: MouseEventInit = {bubbles: true, cancelable: true, composed: true, button: 0,
      clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, ...init};
    // a cancelled pointerdown suppresses the compatibility mouse events, and with mousedown the focus; click fires
    const mouse = el.dispatchEvent(Gesture.pointer('pointerdown', {...base, buttons: 1}));
    if (mouse && el.dispatchEvent(new MouseEvent('mousedown', {...base, buttons: 1})))
      Gesture.focus(el);
    el.dispatchEvent(Gesture.pointer('pointerup', base));
    if (mouse)
      el.dispatchEvent(new MouseEvent('mouseup', base));
    if (init === undefined && el.matches('input[type="checkbox"], input[type="radio"]'))
      el.click();
    else
      el.dispatchEvent(new MouseEvent('click', {...base, detail: 1}));
  }

  static pointer(type: string, init: MouseEventInit): MouseEvent {
    const pointer = {pointerId: 1, pointerType: 'mouse', isPrimary: true, ...init};
    return typeof PointerEvent === 'function' ? new PointerEvent(type, pointer) : new MouseEvent(type, pointer);
  }

  static focus(el: HTMLElement): void {
    const target = el.closest<HTMLElement>(FOCUSABLE);
    if (target !== null)
      target.focus();
    else
      (document.activeElement as HTMLElement | null)?.blur?.();
  }
}

// the field's own selection where it has one; kept here where it has none (a type="number" field, a lenient DOM)
class Caret {
  private constructor(private readonly input: HTMLInputElement, private start: number, private end: number) {}

  static of(input: HTMLInputElement, replace: boolean, focused: boolean): Caret {
    const length = typeof input.value === 'string' ? input.value.length : 0;
    const caret = new Caret(input, replace ? 0 : length, length);
    if (replace || !focused)
      caret._show();
    return caret;
  }

  insert(ch: string): void {
    const input = this.input;
    if (typeof input.value !== 'string') {
      input.append(ch);
      return;
    }
    if (typeof input.selectionStart === 'number' && typeof input.selectionEnd === 'number')
      [this.start, this.end] = [input.selectionStart, input.selectionEnd];
    input.value = input.value.slice(0, this.start) + ch + input.value.slice(this.end);
    this.start = this.end = this.start + 1;
    this._show();
  }

  private _show(): void {
    try {
      this.input.setSelectionRange?.(this.start, this.end);
    }
    catch {
    }
  }
}
