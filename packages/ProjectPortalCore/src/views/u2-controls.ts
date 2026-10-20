/** The views' u2 building blocks: a dialog with named, busy-guarded action buttons, labelled
 * inputs that carry a data-u2-name, a table whose rows can be named, and the discussion composer.
 * u2's Dialog keeps its busy state and error row for OK alone; the package's dialogs have several
 * actions, each of which holds the others while it runs, so the guard stays here. */
import * as grok from 'datagrok-api/grok';
import {BasicTable, BoolInput, ChoiceInput, Control, Dialog, divH, divV, link, signal, span, TextArea,
  TextInput} from '@datagrok-libraries/u2';
import type {MessageInput, ReadonlySignal, TableColumn} from '@datagrok-libraries/u2';
import {messageInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {installStyles} from './styles';

export const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const messageOf = (e: any): string => e?.message ?? String(e);

/** Where a dialog shows what went wrong: a region named [name], each message named [itemName]. */
export class ErrorRegion {
  readonly root: HTMLElement;

  constructor(name: string, private readonly itemName: string) {
    this.root = named(divV([], 'pp-warnings'), name);
  }

  show(e: any): void {
    this.root.replaceChildren(named(span(messageOf(e), 'pp-reason'), this.itemName));
  }

  clear(): void {
    this.root.replaceChildren();
  }
}

export interface DialogButton {
  label: string;
  /** data-u2-name of the button. */
  name: string;
  primary?: boolean;
  /** A throw keeps the dialog open and shows the message; otherwise the dialog closes. */
  run: () => unknown | Promise<unknown>;
  /** True: the dialog stays open after [run] succeeds. */
  keepOpen?: boolean;
  /** After a successful [run] has closed the dialog: what follows, such as a balloon or a redraw. */
  then?: () => unknown | Promise<unknown>;
  enabled?: ReadonlySignal<boolean>;
}

export interface DialogParts {
  content: (HTMLElement | Control)[] | HTMLElement | Control;
  buttons?: DialogButton[];
}

export interface DialogOptions {
  title: string;
  /** data-u2-name of the dialog root. */
  name: string;
  width?: number;
  height?: number;
  /** Where a failed button run is reported; without one, a balloon. */
  errors?: ErrorRegion;
  onClose?: () => void;
}

export type PortalDialog = Dialog & {
  /** Runs the named button as a click would, busy guard included; answers whether it ran and succeeded. */
  press: (name: string) => Promise<boolean>;
};

/** A u2 dialog, shown non-modal as ui.dialog was. [build] runs in the dialog's scope, so every
 * input it makes is disposed with the dialog, which happens once it closes. Action buttons are
 * named, disabled while any one of them runs, and close the dialog when their run succeeds. */
export function openDialog(opts: DialogOptions, build: () => DialogParts): PortalDialog {
  installStyles();
  const dlg = Dialog.create(opts.title, {name: opts.name});
  const busy = signal(false);
  const buttons = new Map<string, HTMLButtonElement>();
  const press = async (name: string): Promise<boolean> => {
    const b = parts.buttons?.find((x) => x.name === name);
    const el = buttons.get(name);
    if (!b || !el || busy.value || el.disabled)
      return false;
    busy.value = true;
    el.classList.add('pp-busy');
    try {
      await b.run();
      opts.errors?.clear();
      if (!b.keepOpen)
        dlg.close();
    } catch (e) {
      if (opts.errors)
        opts.errors.show(e);
      else
        grok.shell.error(messageOf(e));
      return false;
    } finally {
      busy.value = false;
      el.classList.remove('pp-busy');
      // a button disabled while it ran drops the focus to the page, where Escape no longer reaches the dialog
      if (dlg.isOpen.value && !dlg.root.contains(document.activeElement))
        el.focus();
    }
    await Promise.resolve(b.then?.()).catch((e) => grok.shell.error(messageOf(e)));
    return true;
  };
  const parts = dlg.runInScope(build);
  for (const c of Array.isArray(parts.content) ? parts.content : [parts.content])
    dlg.add(c);
  for (const b of parts.buttons ?? []) {
    dlg.addButton(b.label, () => void press(b.name), {primary: b.primary});
    const el = named(dlg.getButton(b.label)!, b.name);
    buttons.set(b.name, el);
    dlg.effect(() => el.disabled = busy.value || (b.enabled ? !b.enabled.value : false));
  }
  // closing disposes the dialog and what it holds
  if (opts.onClose)
    dlg.onClose(opts.onClose);
  dlg.show({width: opts.width, height: opts.height});
  return Object.assign(dlg, {press});
}

/** A single-line text input; [name] lands on the input's root, as the old field wrappers had it. */
export const text = (name: string, label?: string,
  o: {placeholder?: string, value?: string, search?: boolean, inline?: boolean, tooltipText?: string,
    fill?: boolean, onChanged?: (v: string) => void} = {}): TextInput =>
  new TextInput({name, label, ...o});

/** A multi-line text input, [rows] lines tall to start with; it grows with the text. */
export const area = (name: string, label?: string,
  o: {placeholder?: string, value?: string, rows?: number, fill?: boolean} = {}): TextArea =>
  new TextArea({name, label, autoGrow: true, ...o});

/** A native select over {code, name} options. With [blank], an option of that text stands for no
 * choice (value ''); without it, the value starts at [value] or the first option. Read the value
 * with [choiceValue]. */
export function choice(name: string, label: string | undefined, options: {code: string, name: string}[],
  value = '', o: {blank?: string, inline?: boolean, onChanged?: (v: string) => void} = {}): ChoiceInput {
  const items = [...(o.blank !== undefined ? [{value: '', label: o.blank}] : []),
    ...options.map((x) => ({value: x.code, label: x.name}))];
  const start = items.some((i) => i.value === value) ? value : (o.blank !== undefined ? '' : items[0]?.value ?? '');
  return new ChoiceInput({name, label, items, value: start, nullable: false, inline: o.inline,
    onChanged: o.onChanged ? (v) => o.onChanged!(v ?? '') : undefined});
}

/** The chosen code, '' for none. */
export const choiceValue = (input: ChoiceInput): string => input.value.value ?? '';

/** [input] behind a muted [label] on one row, which carries [name] as the old select wrappers did:
 * the filter rows, and the board's mode and persona. */
export const labelled = (name: string, label: string, input: Control): HTMLElement =>
  named(divH([span(label, 'pp-muted'), input.root], 'pp-preference'), name);

/** A checkbox. */
export const check = (name: string, label?: string, value = false,
  o: {enabled?: boolean, tooltipText?: string, onChanged?: (v: boolean) => void} = {}): BoolInput =>
  new BoolInput({name, label, value, ...o});

/** A u2 BasicTable whose rows can be named, classed and titled. Its cells wrap (styles.ts), as the
 * hand-built tables did, so a wide table fits its dialog. */
export function table<T>(name: string, columns: TableColumn<T>[], items: T[],
  rows: {name?: (item: T, i: number) => string | null, cls?: (item: T, i: number) => string | null,
    title?: (item: T, i: number) => string | null} = {}): BasicTable<T> {
  const t = new BasicTable<T>({columns, items, row: (item, i, tr) => {
    const n = rows.name?.(item, i);
    if (n)
      tr.dataset.u2Name = n;
    const c = rows.cls?.(item, i);
    if (c)
      tr.classList.add(...c.split(' '));
    const tip = rows.title?.(item, i);
    if (tip)
      tr.title = tip;
  }});
  t.root.dataset.u2Name = name;
  t.root.classList.add('pp-table');
  return t;
}

/** A compose box: u2's MessageInput with the platform's user mentions. [name] goes on the editable
 * part, `<name>-send` on its Send button. Ctrl+Enter (and Send) call [onSend]; a throw keeps the text.
 * In a dialog, [ownSend] false hides Send: the dialog's button sends instead. */
export function composer(opts: {name: string, placeholder: string, onSend: (markup: string) => Promise<void>,
  ownSend?: boolean, sendLabel?: string}): MessageInput {
  installStyles();
  const input = messageInput({placeholder: opts.placeholder, sendOn: 'ctrlEnter', onSend: opts.onSend});
  input.root.dataset.u2Name = `${opts.name}-composer`;
  const editor = input.root.querySelector<HTMLElement>('.u2-msg-editor');
  if (editor)
    editor.dataset.u2Name = opts.name;
  const send = input.root.querySelector<HTMLElement>('.u2-msg-send');
  if (send) {
    send.dataset.u2Name = `${opts.name}-send`;
    if (opts.sendLabel)
      send.textContent = opts.sendLabel;
  }
  if (opts.ownSend === false)
    input.root.classList.add('pp-msg-nosend');
  return input;
}

/** A u2 link carrying [name]: reachable with Tab, run with a click or Enter. [cls] adds the package's
 * classes (pp-gap ...); [stop] keeps the click from the row it sits in, which opens something else. */
export function linkTo(text: string, name: string, onClick: () => void,
  o: {cls?: string, title?: string, stop?: boolean} = {}): HTMLAnchorElement {
  const el = named(link(text, (e) => {
    if (o.stop)
      e.stopPropagation();
    onClick();
  }), name);
  if (o.cls)
    el.classList.add(...o.cls.split(' '));
  if (o.title)
    el.title = o.title;
  return el;
}
