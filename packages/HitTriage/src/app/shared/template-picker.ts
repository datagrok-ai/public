/* The New Campaign section's Template picker: the app's templates by file name, with New, Clone and Delete beside it,
   and under it the new-campaign form of the template picked — built for that template alone, so a late answer for
   an earlier pick never lands. An app with no template at all offers to create the first one. */
import * as grok from 'datagrok-api/grok';
import {AsyncValue, ChoiceInput, Control, Form, LabelColumn, Scope, batch, button, component, computed, confirmDialog,
  div, each, emptyState, iconButton, signal, skeleton} from '@datagrok-libraries/u2';
import type {AsyncValueState, ReadonlySignal} from '@datagrok-libraries/u2';
import type {TemplateCompute} from '@datagrok-libraries/statistics/src/compute-functions/types';
import {_package} from '../../package';
import {AppName} from '../types';
import {i18n} from '../consts';
import {loadTemplate} from '../utils';

export interface CampaignTemplate {
  name: string;
  key: string;
  compute: TemplateCompute;
}

/** The widest the picker and a template's form get, in px. */
export const PICKER_PX = 720;

/** What a template's form is built with: the label column it shares with the picker, and the scope it lives in (a
 * form built after an await creates its controls there). */
export interface TemplateFormContext {
  column: LabelColumn;
  owner: Scope;
}

export interface TemplatePickerOptions<T extends CampaignTemplate> {
  app: AppName;
  /** A template saved a moment ago: picked first, and offered before the share lists its file. */
  preset?: T;
  /** "+" and Clone: the template editor, empty or filled from a template. */
  edit: (preset?: T) => void;
  form: (template: T, context: TemplateFormContext) => HTMLElement | Control | Promise<HTMLElement | Control>;
}

export interface TemplatePicker {
  /** The templates offered, in the share's order; null until listed. */
  readonly names: ReadonlySignal<string[] | null>;
}

export function templatePicker<T extends CampaignTemplate>(options: TemplatePickerOptions<T>):
  Control & TemplatePicker {
  return component((self) => {
    const names = signal<string[] | null>(null);
    const listing = new AsyncValue(() => list(), {debounceMs: 0});
    self.effect(() => {
      const s = listing.state.value;
      if (s.kind === 'ready')
        names.value = s.value;
    });
    const mode = computed(() => {
      const s = listing.state.value;
      const listed = names.value;
      return s.kind === 'error' ? `error:${s.message}` : listed === null ? 'loading' :
        listed.length === 0 ? 'empty' : 'ready';
    });
    const root = div();
    each(root, () => [mode.value], (m) => m, (m) => section(m));
    return {root, names};

    function section(mode: string): HTMLElement {
      if (mode === 'loading')
        return skeleton(3);
      if (mode.startsWith('error:'))
        return emptyState(`Could not list the templates: ${mode.slice(6)}`, {error: true});
      if (mode === 'empty') {
        return emptyState('No templates yet. A template sets the fields, stages and calculations of every campaign ' +
          'made from it.', {icon: 'file-plus',
          actions: [button('Create template', () => options.edit(), {primary: true})]});
      }
      return picker();
    }

    function picker(): HTMLElement {
      const column = new LabelColumn();
      const form = new Form({labelColumn: column, maxWidth: PICKER_PX});
      // never empty, so no required asterisk
      const input = new ChoiceInput({label: i18n.selectTemplate, name: 'template', nullable: false, required: false,
        items: computed(() => names.value ?? []), value: options.preset?.name ?? names.peek()![0]});
      const picked = new AsyncValue<T>(() => read(input.peek()!),
        {deps: () => input.value, debounceMs: 0});
      input.addOptions(iconButton('plus', () => options.edit(), {tooltip: i18n.createNewTemplate}));
      input.addOptions(iconButton('copy', () => {
        const template = picked.peek();
        if (template !== undefined)
          options.edit(template);
      }, {tooltip: 'Clone template'}));
      input.addOptions(iconButton('trash-alt', () => void remove(input), {tooltip: 'Delete template',
        disabled: computed(() => (names.value?.length ?? 0) < 2 ? 'The last template cannot be deleted' : null)}));
      form.add(input);
      const host = div([]);
      each(host, () => [picked.state.value], (s) => s.kind === 'ready' ? `ready:${s.value.name}` : s.kind,
        (s) => templateForm(s, column));
      return div([form.root, host]);
    }

    function templateForm(state: AsyncValueState<T>, column: LabelColumn): HTMLElement {
      if (state.kind === 'error')
        return emptyState(`Could not read the template: ${state.message}`, {error: true});
      if (state.kind !== 'ready')
        return skeleton(3);
      const owner = Scope.ambient!;
      const built = options.form(state.value, {column, owner});
      if (!(built instanceof Promise))
        return Control.is(built) ? built.root : built;
      const slot = div([skeleton(3)]);
      built.then((el) => {
        if (!owner.isDisposed)
          slot.replaceChildren(Control.is(el) ? el.root : el);
      }, (e) => {
        if (!owner.isDisposed) {
          slot.replaceChildren(emptyState(`Could not build the form: ${e instanceof Error ? e.message : e}`,
            {error: true}));
        }
      });
      return slot;
    }

    async function list(): Promise<string[]> {
      const folder = `${options.app}/templates`;
      const listed = !await _package.files.exists(folder) ? [] : (await _package.files.list(folder))
        .filter((file) => file.name.endsWith('.json'))
        .map((file) => file.name.slice(0, -5));
      const preset = options.preset?.name;
      return preset !== undefined && !listed.includes(preset) ? [...listed, preset] : listed;
    }

    async function read(name: string): Promise<T> {
      const preset = options.preset;
      return preset?.name === name ? preset : await loadTemplate<T>(`${options.app}/templates/${name}.json`);
    }

    async function remove(input: ChoiceInput): Promise<void> {
      const name = input.peek();
      if (name === null)
        return;
      const deleted = await confirmDialog({title: 'Delete template', name: 'delete-template', confirmText: 'Delete',
        destructive: true, message: `Are you sure you want to delete template ${name}?`,
        detail: 'Campaigns made from it keep working.',
        onConfirm: () => _package.files.delete(`${options.app}/templates/${name}.json`)});
      if (!deleted)
        return;
      const rest = (names.peek() ?? []).filter((n) => n !== name);
      batch(() => {
        names.value = rest;
        input.value = rest[0] ?? null;
      });
      grok.shell.info(`Template ${name} was deleted successfully`);
    }
  }, {name: 'template-picker'});
}
