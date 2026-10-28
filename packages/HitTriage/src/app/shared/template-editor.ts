/* The New Template view of Hit Design and PeptiHit: one page, in this order — Details (Name, Key, Layout),
   Additional fields, Stages, Compute, Submit — over a footer that stays in view, with a line saying what CREATE makes
   or why it cannot yet, CANCEL and CREATE. A clone starts from a template under a name and a key no template has, and
   keeps its fields, stages, compute, submit function and layout. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BoolInput, ChoiceInput, Control, Form, Page, RecordsInput, Section, TextInput, batch, button, component,
  computed, confirmDialog, div, on, signal, statusLine, untracked} from '@datagrok-libraries/u2';
import {ViewHost, fromDartInput} from '@datagrok-libraries/u2/dg';
import {_package} from '../../package';
import {AppName, CampaignFieldTypes, ComputeFunctions, HitDesignTemplate, HitTriageCampaignField,
  TemplateCompute} from '../types';
import {HitTriageSubmitTag, HTFunctionOrderingLSKey, i18n} from '../consts';
import {getFuncPackageNameSafe} from '../utils';
import {descriptorsTree, rerunOnOpenOption} from '../dialogs/functions-dialog';
import {TemplateLayout, layoutFile, layoutFileInput, readLayout} from '../accordeons/layout-input';
import {computeEditor} from './compute-editor';
import type {HitAppBase} from '../hit-app-base';

export interface TemplateEditorRequest<T extends HitDesignTemplate> {
  app: HitAppBase<unknown>;
  /** Clone: the template the new one starts from. */
  preset?: T;
  /** Runs once the template file is written, before the editor closes. */
  onCreated: (template: T) => unknown;
}

/** Opens the New Template view over the current one, the app's landing. Its crumb leads back there and closes the
 * editor, as it always did — after asking, once something was typed into it. */
export async function openTemplateEditor<T extends HitDesignTemplate>(request: TemplateEditorRequest<T>):
  Promise<void> {
  const {app, preset} = request;
  const parent = grok.shell.v;
  ui.setUpdateIndicator(parent.root, true);
  try {
    const existing = await existingTemplates(app.appName);
    let closing = false;
    const close = () => {
      closing = true;
      view.close();
    };
    let editor!: Control & TemplateEditor;
    const view = ViewHost.create({name: i18n.createNewTemplate, content: () => {
      editor = templateEditor<T>({app: app.appName, preset, existing,
        functions: app.computeFunctions, cancel: close,
        created: async (template) => {
          try {
            await request.onCreated(template);
          }
          catch (e) {
            grok.shell.error(`Template ${template.name} was created, but the app could not show it: ` +
              `${e instanceof Error ? e.message : e}`);
          }
          close();
        }});
      return {root: editor, path: `${new URL(app.baseUrl).pathname}/new-template`};
    }}).view;
    view.parentView = parent;
    view.parentCall = app.parentCall;
    view.helpUrl = `${_package.webRoot}README_HD.md#templates`;
    grok.shell.addView(view);
    const leave = async () => {
      if (!editor.dirty || await confirmDialog({title: 'Discard the new template?', name: 'discard-template',
        message: 'What you filled in for the new template is not saved.', confirmText: 'Discard', destructive: true}))
        close();
      else
        grok.shell.v = view;
    };
    // the crumb (or any way back) to the landing leaves the editor; another view leaves it open behind
    let current = true;
    editor.subscribe(grok.events.onCurrentViewChanged, () => {
      const now = grok.shell.v;
      const left = current && !closing && now?.dart === parent.dart &&
        Array.from(grok.shell.views).some((v) => v.dart === view.dart);
      current = now?.dart === view.dart;
      if (left)
        void leave();
    });
  }
  catch (e) {
    grok.shell.error(`Could not open the template editor: ${e instanceof Error ? e.message : e}`);
  } finally {
    ui.setUpdateIndicator(parent.root, false);
  }
}

interface TemplateId {
  name: string;
  key: string;
}

export interface TemplateEditorOptions<T extends HitDesignTemplate> {
  app: AppName;
  preset?: T;
  /** The app's templates: a new one takes a name and a key none of them has. */
  existing: TemplateId[];
  functions: ComputeFunctions;
  cancel: () => void;
  created: (template: T) => Promise<unknown>;
}

const NO_COMPUTE: TemplateCompute = {descriptors: {enabled: false, args: []}, functions: [], scripts: [], queries: []};
const FIELD_TYPES = Object.keys(CampaignFieldTypes);

export interface TemplateEditor {
  /** Something was typed or picked: an input edited, a row added, removed or moved, a descriptor picked. */
  readonly dirty: boolean;
}

export function templateEditor<T extends HitDesignTemplate>(options: TemplateEditorOptions<T>):
  Control & TemplateEditor {
  return component((self) => {
    // the names and keys taken, trimmed and in lower case; a key by the template that has it
    const names = new Set<string>();
    const keys = new Map<string, string>();
    const submitFuncs = submitFunctions();
    const reading = signal(false);
    const busy = signal(false);
    const edited = signal(false);
    let reads = 0;
    remember(options.existing);
    const preset = options.preset;
    const layout = signal<TemplateLayout | null>(preset?.layoutViewState ?
      {viewState: preset.layoutViewState, localPath: preset.localLayoutPath ?? null} : null);

    const form = new Form();
    const nameInput = new TextInput({label: 'Name', name: 'name', required: true, tooltipText: 'Template name',
      value: preset ? freeName(preset.name, names) : ''});
    nameInput.addValidator((v) => names.has(v.trim().toLowerCase()) ?
      'A template with this name already exists' : null);
    const keyInput = new TextInput({label: 'Key', name: 'key', required: true,
      tooltipText: 'Template key used for campaign prefix',
      value: preset ? freeKey(preset.key, new Set(keys.keys())) : ''});
    keyInput.addValidator((v) => {
      const owner = keys.get(v.trim().toLowerCase());
      return owner === undefined ? null : `${v.trim()} is used by ${owner}`;
    });

    const picker = layoutFileInput();
    picker.setTooltip('Select a layout file');
    const carried = layout.peek();
    if (carried !== null)
      picker.value = layoutFile(carried, preset!.name);
    const layoutInput = fromDartInput(picker, 'layout', {fill: true});
    layoutInput.addValidator((file) => file !== null && !reading.value && layout.value === null ?
      'Not a layout file' : null);
    layoutInput.effect(() => {
      const file = layoutInput.value;
      untracked(() => void readLayoutFile(file));
    });

    const fields = new RecordsInput<HitTriageCampaignField>({label: '', name: 'fields', addText: 'Add field',
      value: preset?.campaignFields ?? [], newItem: () => ({name: '', type: 'String', required: false}),
      rowName: (f, i) => f.name || `Field ${i + 1}`,
      columns: [
        {name: 'name', header: 'Name', input: (v) => new TextInput({inline: true, value: v ?? ''})},
        {name: 'type', header: 'Type', width: '120px',
          input: (v) => new ChoiceInput({inline: true, nullable: false, items: FIELD_TYPES, value: v ?? 'String'})},
        {name: 'required', header: 'Required', width: '64px',
          input: (v) => new BoolInput({inline: true, value: v === true})},
      ]});
    const stages = new RecordsInput<{name: string}>({label: '', name: 'stages', addText: 'Add stage',
      orderable: true, value: (preset?.stages ?? []).map((name) => ({name})), newItem: () => ({name: ''}),
      rowName: (_, i) => `Stage ${i + 1}`,
      columns: [{name: 'name', header: 'Name', input: (v) => new TextInput({inline: true, value: v ?? ''})}]});

    const submit = preset?.submit;
    const presetSubmit = submit === undefined ? undefined : submitFuncs.find((f) => f.name === submit.fName &&
      (!submit.package || getFuncPackageNameSafe(f) === submit.package));
    const submitInput = new ChoiceInput({label: 'Submit function', name: 'submit-function', nullable: true,
      tooltipText: 'Select function to be called upon submitting', value: presetSubmit?.nqName ?? null,
      items: submitFuncs.map((f) => ({value: f.nqName, label: f.friendlyName ?? f.name}))});

    const compute = computeEditor({functions: options.functions, value: preset?.compute ?? NO_COMPUTE,
      descriptors: descriptorsTree, orderKey: HTFunctionOrderingLSKey, rerunOnOpen: rerunOnOpenOption});

    const section = (title: string, hint?: string): HTMLElement => {
      const s = new Section({title, hint, collapsible: false});
      form.addElement(s.root);
      return s.body;
    };
    const details = section('Details');
    form.add(nameInput, details).add(keyInput, details).add(layoutInput, details);
    form.add(fields, section('Additional fields', 'filled in by people when they start a campaign'));
    form.add(stages, section('Stages', 'new molecules start in the first stage'));
    section('Compute', 'calculated for every molecule').append(compute.root);
    form.add(submitInput, section('Submit'));

    const missing = computed(() => [nameInput, keyInput].filter((i) => i.value.trim() === '')
      .map((i) => i.label));
    const quiet = computed(() => {
      const key = keyInput.value.trim();
      return missing.value.length > 0 ? `Fill in ${missing.value.join(' and ')} to create the template` :
        `Creates template "${nameInput.value.trim()}"; campaigns get codes ${key}-1, ${key}-2, …`;
    });
    const create = button(i18n.createTemplate, () => void createTemplate(), {primary: true,
      disabled: computed(() => busy.value ? 'Creating the template…' :
        reading.value ? 'Reading the layout…' : form.problem.value)});
    const page = new Page({content: form, maxWidth: 820, labelColumn: false});
    form.addButtons((row) => row.append(statusLine(form.shownProblem, quiet),
      button(i18n.cancel, () => options.cancel()), create), {host: page.footer});
    const root = div([page], 'hit-triage-template-editor');

    // typing and picking reach the root as input and change events; rows added, removed or moved do not
    for (const type of ['input', 'change'])
      on(root, type, () => edited.value = true);
    let first = true;
    self.effect(() => {
      fields.value;
      stages.value;
      if (!first)
        untracked(() => edited.value = true);
      first = false;
    });
    return {root, get dirty() {
      return edited.peek() || compute.edited.peek();
    }};

    /** The template the editor holds, its keys in the order the editor has always written them. */
    function currentTemplate(): T {
      const chosen = layout.peek();
      const submitFunc = submitFuncs.find((f) => f.nqName === submitInput.peek());
      return {
        name: nameInput.peek().trim(),
        key: keyInput.peek().trim(),
        campaignFields: fields.peek().filter((f) => f.name).map((f) => ({...f, required: !!f.required})),
        stages: stages.peek().map((s) => s.name).filter((s) => s),
        ...chosen === null ? {} :
          {layoutViewState: chosen.viewState, ...chosen.localPath === null ? {} : {localLayoutPath: chosen.localPath}},
        compute: compute.read(),
        ...submitFunc === undefined ? {} :
          {submit: {fName: submitFunc.name, package: getFuncPackageNameSafe(submitFunc)}},
      } as unknown as T;
    }

    async function readLayoutFile(file: DG.FileInfo | null): Promise<void> {
      const read = ++reads;
      reading.value = true;
      const loaded = file === null ? null : await readLayout(file);
      if (read !== reads || self.isDisposed)
        return;
      batch(() => {
        layout.value = loaded;
        reading.value = false;
      });
      layoutInput.revalidate();
    }

    function remember(existing: TemplateId[]): void {
      for (const t of existing) {
        if (t.name.trim() !== '')
          names.add(t.name.trim().toLowerCase());
        if (t.key.trim() !== '')
          keys.set(t.key.trim().toLowerCase(), t.name);
      }
    }

    /** Writes the template after reading the app's templates again: a name or a key another session took since the
     * editor opened, or a failed write, keeps the editor open. CREATE stays busy after a write: the editor closes. */
    async function createTemplate(): Promise<void> {
      const template = currentTemplate();
      busy.value = true;
      try {
        remember(await existingTemplates(options.app));
        const taken = [nameInput, keyInput].filter((input) => {
          input.revalidate();
          return input.validity.peek() !== null;
        });
        if (taken.length > 0) {
          busy.value = false;
          for (const input of taken)
            input.touch();
          return;
        }
        await _package.files.writeAsText(`${options.app}/templates/${template.name}.json`, JSON.stringify(template));
      }
      catch (e) {
        busy.value = false;
        grok.shell.error(`Could not save template ${template.name}: ${e instanceof Error ? e.message : e}`);
        return;
      }
      grok.shell.info('Template created successfully');
      await options.created(template);
    }
  }, {name: 'template-editor'});
}

/** The app's templates by name and key; one that cannot be read, by name alone. */
async function existingTemplates(app: AppName): Promise<TemplateId[]> {
  const folder = `${app}/templates`;
  if (!await _package.files.exists(folder))
    return [];
  const files = (await _package.files.list(folder)).filter((f) => f.name.endsWith('.json'));
  return Promise.all(files.map(async (f) => {
    const name = f.name.slice(0, -5);
    try {
      return {name, key: String(JSON.parse(await _package.files.readAsText(f)).key ?? '')};
    }
    catch {
      return {name, key: ''};
    }
  }));
}

/** The functions a campaign can be submitted to, each once. */
function submitFunctions(): DG.Func[] {
  const found = DG.Func.find({meta: {role: HitTriageSubmitTag}}).concat(DG.Func.find({tags: [HitTriageSubmitTag]}));
  return Array.from(new Map(found.map((f) => [f.nqName, f])).values());
}

/** "<name> (copy)", or "(copy 2)" and on when that is taken too; `taken` in lower case. */
export function freeName(name: string, taken: Set<string>): string {
  for (let n = 1; ; n++) {
    const candidate = n === 1 ? `${name} (copy)` : `${name} (copy ${n})`;
    if (!taken.has(candidate.toLowerCase()))
      return candidate;
  }
}

/** The key with the next free number, as wide as its own: KHL → KHL2, KHL2 → KHL3, KHL09 → KHL10; none for none. */
export function freeKey(key: string, taken: Set<string>): string {
  if (key === '')
    return '';
  const [, base, digits] = /^(.*?)(\d*)$/.exec(key)!;
  const next = (n: number) => `${base}${String(n).padStart(digits.length, '0')}`;
  let n = digits === '' ? 2 : Number(digits) + 1;
  while (taken.has(next(n).toLowerCase()))
    n++;
  return next(n);
}
