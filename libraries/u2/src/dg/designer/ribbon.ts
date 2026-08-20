/* The ribbon and its dialogs: New, the Open menu (samples, gallery, paste), Copy, undo/redo, the
   editing verbs and the bounds toggle — and every finished action hands focus back to the canvas. */
import * as grok from 'datagrok-api/grok';
import {Scope} from '../../core/scope.js';
import type {ReadonlySignal, Signal} from '../../core/signals.js';
import type {Action} from '../../components/actions/actions.js';
import {dropDownButton, iconButton} from '../../components/actions/buttons.js';
import {Dialog} from '../../components/containers/dialog.js';
import type {Menu} from '../../components/navigation/menu.js';
import {TextArea, TextInput} from '../../components/inputs/text-input.js';
import type {SpecEditor} from '../../spec/editor.js';
import {ACTIONS, ACTION_ICONS} from './actions.js';
import {loadGallery, saveToGallery} from './gallery.js';
import {SAMPLES} from './samples.js';

export interface RibbonHost {
  readonly mode: HTMLElement;
  readonly outlines: Signal<boolean>;
  readonly editor: ReadonlySignal<SpecEditor | undefined>;
  effect(fn: () => void): void;
  revision(): number;
  actions(): Action[];
  run(name: string): void;
  open(spec: string | object): void;
  dump(): object | null;
  report(): void;
  refocus(): void;
}

export class Ribbon {
  constructor(private readonly _host: RibbonHost) {}

  build(): HTMLElement[] {
    const host = this._host;
    const undo = iconButton('undo', () => {
      host.editor.peek()?.undo();
      host.refocus();
    }, {tooltip: 'Undo'});
    const redo = iconButton('redo', () => {
      host.editor.peek()?.redo();
      host.refocus();
    }, {tooltip: 'Redo'});
    host.effect(() => undo.disabled = host.editor.value?.canUndo.value !== true);
    host.effect(() => redo.disabled = host.editor.value?.canRedo.value !== true);
    // `actions()` answers these verbs in this order, or nothing at all when there is no
    // selection to act on — so the buttons read their enablement back by position
    const actions = ACTIONS.map((name) =>
      iconButton(ACTION_ICONS[name], () => host.run(name), {tooltip: name}));
    host.effect(() => {
      const current = host.revision() >= 0 ? host.actions() : [];
      for (let i = 0; i < actions.length; i++)
        actions[i].disabled = current[i] === undefined || current[i].enabled === false;
    });
    return [
      // first, not last: the platform ribbon never reflows, and a narrow view container clips
      // the tail — the casualties must be icon verbs (all in the menus too), never the mode
      host.mode,
      iconButton('file', () => this.newForm(), {tooltip: 'New form'}),
      // one button, one menu: samples, gallery and the paste dialog all live under Open —
      // a split button hid the samples behind an unnamed chevron
      dropDownButton('', (menu) => this._openMenu(menu), {icon: 'folder-open', tooltip: 'Open'}),
      iconButton('copy', () => {
        this._copy();
        host.refocus();
      }, {tooltip: 'Copy spec'}),
      undo,
      redo,
      ...actions,
      iconButton('border-all', () => host.refocus(),
        {toggle: host.outlines, tooltip: 'Control bounds'}),
    ];
  }

  /** F3/F12 loads share one gate: `dirty` is what the status bar's `modified` reads — undoing back
   * to the last save is clean, and a saved document keeps its history. */
  private _load(spec: object): void {
    const proceed = () => {
      this._host.open(spec);
      this._host.refocus();
    };
    if (this._host.editor.peek()?.dirty.peek() !== true) {
      proceed();
      return;
    }
    const scope = new Scope();
    Scope.runWith(scope, () => new Dialog('Discard changes?'))
      .add('The form has been modified. Discard the changes?')
      .onOK(() => {
        scope.dispose();
        proceed();
      })
      .onCancel(() => {
        scope.dispose();
        this._host.refocus();
      })
      .show({modal: true});
  }

  newForm(): void {
    this._load(SAMPLES[0].spec);
  }

  /** The Open arrow's menu, rebuilt per open so the gallery names never show stale state. */
  private _openMenu(menu: Menu): void {
    menu.group('Samples');
    for (const sample of SAMPLES)
      menu.item(sample.name, () => this._load(sample.spec));
    menu.endGroup();
    const saved = loadGallery();
    const names = Object.keys(saved);
    if (names.length > 0) {
      menu.group('Gallery');
      for (const name of names)
        menu.item(name, () => this._load(saved[name]));
      menu.endGroup();
    }
    menu.separator();
    menu.item('Open spec…', () => this._openDialog());
    menu.item('Save to gallery…', () => this.saveToGallery());
  }

  saveToGallery(): void {
    const spec = this._host.dump();
    if (!spec) {
      grok.shell.info('Nothing to save — the canvas is empty.');
      this._host.refocus();
      return;
    }
    const scope = new Scope();
    const name = Scope.runWith(scope, () => new TextInput({label: 'Name', placeholder: 'My form'}));
    Scope.runWith(scope, () => new Dialog('Save to gallery'))
      .add(name)
      .add('An entry that already carries this name is overwritten.')
      .onOK(() => {
        const trimmed = name.value.peek().trim();
        scope.dispose();
        if (trimmed === '') {
          grok.shell.warning('A gallery entry needs a name.');
          this._host.refocus();
          return;
        }
        if (saveToGallery(trimmed, spec)) {
          grok.shell.info(`Saved to the gallery as "${trimmed}".`);
          this._host.editor.peek()?.markClean();
          this._host.report();
        } else
          grok.shell.warning(`Could not save "${trimmed}" — storage is unavailable or full.`);
        this._host.refocus();
      })
      .onCancel(() => {
        scope.dispose();
        this._host.refocus();
      })
      .show({modal: true, width: 360});
  }

  private _openDialog(): void {
    const scope = new Scope();
    const spec = this._host.dump();
    const area = Scope.runWith(scope, () => new TextArea({
      label: 'Spec',
      value: spec ? JSON.stringify(spec, null, 2) : '',
      placeholder: 'Paste a dg-ui/1 spec',
    }));
    area.root.classList.add('u2-designer-spec-editor');
    Scope.runWith(scope, () => new Dialog('Open spec'))
      .add(area)
      .onOK(() => {
        this._host.open(area.value.peek());
        scope.dispose();
        this._host.refocus();
      })
      .onCancel(() => {
        scope.dispose();
        this._host.refocus();
      })
      .show({modal: true, width: 560});
  }

  private _copy(): void {
    const spec = this._host.dump();
    if (!spec) {
      grok.shell.info('Nothing to copy — the canvas is empty.');
      return;
    }
    navigator.clipboard.writeText(JSON.stringify(spec, null, 2))
      .then(() => grok.shell.info('Spec copied to the clipboard.'),
        (e: unknown) => grok.shell.error(`Could not copy: ${e instanceof Error ? e.message : String(e)}`));
  }
}
