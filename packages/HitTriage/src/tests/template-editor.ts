/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import {loadTemplate} from '../app/utils';
import {openTemplateEditor} from '../app/shared/template-editor';
import {createTemplateAccordeon} from '../app/accordeons/new-template-accordeon';
import {createPepTriageTemplateAccordeon} from '../app/accordeons/new-pep-triage-template-accordeon';
import {HitAppBase} from '../app/hit-app-base';
import {HitDesignTemplate, IComputeDialogResult, TemplateCompute} from '../app/types';
import {ComputeEditorOptions, computeEditor, templateCompute} from '../app/shared/compute-editor';
import {chemFunctionsDialog, descriptorsTree, rerunOnOpenOption} from '../app/dialogs/functions-dialog';
import {HTFunctionOrderingLSKey} from '../app/consts';
import {expectSame, noComputeTemplate, pepTriageTemplate, Sandbox, triageTemplate} from './fixtures';
import {onRoot} from './owners';

/** A name the server refuses as a file name (a NUL in it: 400, on any storage): the write fails, as a full disk or a
 * lost connection would make it. */
const unwritable = (app: string) => `ZZB3 ${app} \u0000`;
const saveFailed = (app: string) => `Could not save template ZZB3 ${app}`;

const openViews = () => Array.from(grok.shell.views).map((v) => v.name);

const computeEditorOnRoot = (options: ComputeEditorOptions) => onRoot('compute editor test', () => computeEditor(options));

/** A compute function, script and query found by name; a stand without one fails the test rather than checks less. */
function required<T>(found: T | undefined, what: string): T {
  if (found === undefined)
    throw new Error(`the stand has no ${what}`);
  return found;
}

interface OpenEditor<T> {
  root: HTMLElement;
  /** What the editor handed to the landing after writing it. */
  created: Promise<T>;
  /** Whether it did. */
  isCreated: () => boolean;
}

/** The Hit Design editor opened as the template picker opens it, once its Compute section has its descriptors. */
async function openEditor<T extends HitDesignTemplate>(app: HitAppBase<unknown>, preset?: T): Promise<OpenEditor<T>> {
  let created!: (t: T) => void;
  let isCreated = false;
  const done = new Promise<T>((resolve) => created = resolve);
  await openTemplateEditor<T>({app, preset, onCreated: (t) => {
    isCreated = true;
    created(t);
  }});
  await waitUntil(() => grok.shell.v?.name === 'New Template', 'the template editor did not open', 15000);
  const root = grok.shell.v.root.querySelector<HTMLElement>('[data-u2-name="template-editor"]')!;
  await waitUntil(() => root.querySelector('.hit-triage-compute-descriptors .u2-tree:not([hidden])') !== null,
    'the Compute section never listed the descriptors', 30000);
  return {root, created: done, isCreated: () => isCreated};
}

function typeInto(input: HTMLInputElement, text: string): void {
  input.value = text;
  input.dispatchEvent(new Event('input', {bubbles: true}));
}

const editorInput = (root: HTMLElement, name: string) => root.querySelector<HTMLInputElement>(`[data-u2-name="${name}"] input`)!;
const createButton = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLButtonElement>('button'))
  .find((b) => b.textContent?.trim().toLowerCase() === 'create')!;
const disabled = (button: HTMLButtonElement) => button.getAttribute('aria-disabled') === 'true' || button.classList.contains('d4-disabled');
const balloonShown = (text: string) => Array.from(document.querySelectorAll('.d4-balloon')).some((b) => b.textContent?.includes(text));

async function closeEditor(): Promise<void> {
  if (grok.shell.v?.name === 'New Template')
    grok.shell.v.close();
}

/** An input of the Hit Triage / PepTriage editor by its caption, set as typing sets it. */
function typeTriage(root: HTMLElement, caption: string, text: string): void {
  typeInto(root.querySelector<HTMLInputElement>(`[name="input-host-${caption}"] input, ` +
    `[name="input-host-${caption.replace(/ /g, '-')}"] input`)!, text);
}

/** The Hit Triage editor reads and writes the app's own folder, not a sandbox's: what a test puts there goes after it. */
const TRIAGE_TEMPLATES = 'Hit Triage/templates';

async function removeTriageTemplate(name: string): Promise<void> {
  const path = `${TRIAGE_TEMPLATES}/${name}.json`;
  if (await _package.files.exists(path).catch(() => false))
    await _package.files.delete(path);
}

category('Template editor', () => {
  test('Hit Design: a clone proposes a free name and key and keeps the fields, stages, compute, submit function and layout', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const viewState = DG.ViewLayout.fromJson(await _package.files.readAsText('Hit Design/HD-template-layout.layout')).viewState;
      const script = required(DG.Func.find({package: 'ApiSamples', name: 'Pearson'})[0] as DG.Script | undefined,
        'script ApiSamples:Pearson (a table and a column, as compute scripts take them)');
      const query = required(DG.Func.find({package: 'Chembl', name: 'NamesToSmiles'})[0] as DG.DataQuery | undefined,
        'query Chembl:NamesToSmiles (a list first, as compute queries take it)');
      const t: HitDesignTemplate = {...noComputeTemplate('ZZB3 Source', 'ZZB3'), layoutViewState: viewState,
        compute: {descriptors: {enabled: true, args: ['MolWt']}, functions: [],
          scripts: [{name: script.name, id: script.id, args: {}, rerunOnOpen: true}],
          queries: [{name: query.name, id: query.id, args: {}}]}};
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      // the stand tags no script or query as a compute function: the app is told of these two
      app.computeFunctions.scripts.push(script);
      app.computeFunctions.queries.push(query);
      const editor = await openEditor(app, t);
      expect(editorInput(editor.root, 'name').value, 'ZZB3 Source (copy)', 'the proposed name');
      expect(editorInput(editor.root, 'key').value, 'ZZB4', 'the proposed key');
      expect(editor.root.querySelector<HTMLSelectElement>('[data-u2-name="submit-function"] select')!.selectedOptions[0]?.textContent,
        'Demo File Submit', 'the submit function, by its friendly name');
      expect(editorInput(editor.root, 'layout').value, 'ZZB3 Source.layout', 'the layout, shown by its file name');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed disabled for the clone', 10000);
      createButton(editor.root).click();
      const saved = await editor.created;
      const expected = {...t, name: 'ZZB3 Source (copy)', key: 'ZZB4',
        campaignFields: t.campaignFields.map((f) => ({...f, required: !!f.required}))};
      // the arguments of a script or a query are what their forms hold, every input past the molecules
      const withoutArgs = (x: HitDesignTemplate) => ({...x, compute: {...x.compute,
        scripts: x.compute.scripts?.map((s) => ({...s, args: null})),
        queries: x.compute.queries?.map((q) => ({...q, args: null}))}});
      const expectedCompute = {...expected.compute, scripts: [{...t.compute.scripts![0]}], queries: [{...t.compute.queries![0], rerunOnOpen: false}]};
      expectSame(withoutArgs(saved), withoutArgs({...expected, compute: expectedCompute}), 'the clone, its script and query included');
      expectSame(withoutArgs(await loadTemplate(`${sb.appName}/templates/ZZB3 Source (copy).json`)),
        withoutArgs({...expected, compute: expectedCompute}), 'the clone on the share');
      expectSame(Object.keys(saved), ['name', 'key', 'campaignFields', 'stages', 'layoutViewState', 'compute', 'submit'],
        'the keys, in the order the editor has always written them');
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a clone of a template whose layout is a file on a share keeps that file\'s path', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const json = await _package.files.readAsText('Hit Design/HD-template-layout.layout');
      await _package.files.writeAsText(`${sb.appName}/layouts/zz.layout`, json);
      const share = `System:AppData/HitTriage/${sb.appName}/layouts/zz.layout`;
      const t: HitDesignTemplate = {...noComputeTemplate('ZZB3 Shared Layout', 'ZZBS'),
        layoutViewState: DG.ViewLayout.fromJson(json).viewState, localLayoutPath: share};
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      const editor = await openEditor(app, t);
      expect(editorInput(editor.root, 'layout').value.endsWith('zz.layout'), true, 'the layout, shown by its file');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed disabled for the clone', 10000);
      createButton(editor.root).click();
      const saved = await editor.created;
      expect(saved.localLayoutPath, share, 'the share path, re-read by "Apply template layout"');
      expect(saved.layoutViewState, t.layoutViewState, 'the layout');
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: the crumb back to the landing closes the editor, and asks first once something was typed', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const discard = () => document.querySelector<HTMLElement>('[data-u2-name="discard-template"]');
      const press = (caption: string) => Array.from(discard()!.querySelectorAll<HTMLButtonElement>('button'))
        .find((b) => b.textContent?.trim().toLowerCase() === caption)!.click();
      await openEditor(app);
      grok.shell.v = app.multiView;
      await waitUntil(() => !openViews().includes('New Template'), 'the untouched editor stayed open', 10000);
      expect(discard(), null, 'a question about an untouched editor');

      const editor = await openEditor(app);
      typeInto(editorInput(editor.root, 'name'), 'ZZB3 Typed');
      grok.shell.v = app.multiView;
      await waitUntil(() => discard() !== null, 'nothing asked before the typed name went', 10000);
      press('cancel');
      await waitUntil(() => grok.shell.v?.name === 'New Template', 'CANCEL did not bring the editor back', 10000);
      expect(editorInput(editor.root, 'name').value, 'ZZB3 Typed', 'the typed name, kept');
      grok.shell.v = app.multiView;
      await waitUntil(() => discard() !== null, 'nothing asked the second time', 10000);
      press('discard');
      await waitUntil(() => !openViews().includes('New Template'), 'DISCARD left the editor open', 10000);
      expect(grok.shell.v?.name, app.multiView.name, 'the landing, current');
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: with no template at all, the landing\'s editor creates the first one and offers it', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const picker = () => app.infoView.root.querySelector<HTMLElement>('[data-u2-name="template-picker"]');
      const offer = Array.from(picker()?.querySelectorAll<HTMLButtonElement>('[data-u2="empty-state"] button') ?? [])
        .find((b) => b.textContent?.trim() === 'Create template');
      expect(offer !== undefined, true, 'the picker\'s Create template');
      offer!.click();
      await waitUntil(() => grok.shell.v?.name === 'New Template', 'the template editor did not open', 15000);
      const root = grok.shell.v.root.querySelector<HTMLElement>('[data-u2-name="template-editor"]')!;
      typeInto(editorInput(root, 'name'), 'ZZB3 First');
      typeInto(editorInput(root, 'key'), 'ZZBF');
      await waitUntil(() => !disabled(createButton(root)), 'CREATE stayed disabled', 30000);
      createButton(root).click();
      await waitUntil(() => grok.shell.v?.name !== 'New Template' &&
        picker()?.querySelector<HTMLSelectElement>('select')?.value === 'ZZB3 First', 'the landing never offered the first template', 15000);
      expect((await loadTemplate<HitDesignTemplate>(`${sb.appName}/templates/ZZB3 First.json`)).key, 'ZZBF', 'the key on the share');
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: CREATE is disabled while the name or the key is empty or taken', async () => {
    const sb = await Sandbox.create('editor');
    try {
      await sb.writeTemplate(noComputeTemplate('ZZB3 Taken', 'ZZBT'));
      const app = await sb.hitDesign();
      const {root} = await openEditor(app);
      expect(disabled(createButton(root)), true, 'CREATE with nothing filled in');
      typeInto(editorInput(root, 'name'), '   ');
      typeInto(editorInput(root, 'key'), 'ZZBF');
      await waitUntil(() => disabled(createButton(root)), 'CREATE enabled for a name of spaces', 5000);
      typeInto(editorInput(root, 'name'), 'ZZB3 Free');
      typeInto(editorInput(root, 'key'), '  ');
      await waitUntil(() => disabled(createButton(root)), 'CREATE enabled for a key of spaces', 5000);
      typeInto(editorInput(root, 'key'), 'ZZBF');
      typeInto(editorInput(root, 'name'), 'zzb3 taken ');
      await waitUntil(() => disabled(createButton(root)), 'CREATE enabled for a taken name (another case, a space after)', 5000);
      typeInto(editorInput(root, 'name'), 'ZZB3 Free');
      await waitUntil(() => !disabled(createButton(root)), 'CREATE stayed disabled with a free name and key', 5000);
      typeInto(editorInput(root, 'key'), 'zzbt');
      await waitUntil(() => disabled(createButton(root)), 'CREATE enabled for a taken key', 5000);
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a name or a key another session took meanwhile is refused, nothing overwritten or duplicated', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const editor = await openEditor(app);
      typeInto(editorInput(editor.root, 'name'), 'ZZB3 Late');
      typeInto(editorInput(editor.root, 'key'), 'ZZBL');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed disabled', 5000);
      const other = noComputeTemplate('ZZB3 Late', 'ZZBO');
      await sb.writeTemplate(other);
      createButton(editor.root).click();
      await waitUntil(() => disabled(createButton(editor.root)) && editor.root.querySelector('[data-u2-name="name"] [data-u2-part="error"]')
        ?.textContent === 'A template with this name already exists', 'the name was not refused', 10000);
      expect(editor.isCreated(), false, 'the editor reported a template created');
      expect(grok.shell.v.name, 'New Template', 'the view after CREATE');
      expectSame(await loadTemplate(`${sb.appName}/templates/ZZB3 Late.json`), other, 'the other session\'s template');

      typeInto(editorInput(editor.root, 'name'), 'ZZB3 Later');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed disabled with a free name', 5000);
      await sb.writeTemplate(noComputeTemplate('ZZB3 Other', 'ZZBL'));
      createButton(editor.root).click();
      await waitUntil(() => disabled(createButton(editor.root)) && editor.root.querySelector('[data-u2-name="key"] [data-u2-part="error"]')
        ?.textContent === 'ZZBL is used by ZZB3 Other', 'the key was not refused', 10000);
      expect(editor.isCreated(), false, 'the editor reported a template created');
      expect(await _package.files.exists(`${sb.appName}/templates/ZZB3 Later.json`), false, 'a second template with the key');
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a template the share refuses keeps the editor open and says why', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const editor = await openEditor(app);
      typeInto(editorInput(editor.root, 'name'), unwritable('design'));
      typeInto(editorInput(editor.root, 'key'), 'ZZBU');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed disabled', 5000);
      DG.Balloon.closeAll();
      createButton(editor.root).click();
      await waitUntil(() => balloonShown(saveFailed('design')), 'no balloon said the save failed', 15000);
      expect(editor.isCreated(), false, 'the editor reported a template created');
      expect(balloonShown('Template created successfully'), false, 'a success balloon');
      expect(grok.shell.v.name, 'New Template', 'the view after the failed save');
      await waitUntil(() => !disabled(createButton(editor.root)), 'CREATE stayed busy after the failure', 5000);
    } finally {
      await closeEditor();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Compute: the editor reads a template\'s compute back as the compute dialog does', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const toxicity = required(app.computeFunctions.functions.find((f) => f.name === 'addChemRisksColumns'),
        'compute function Chem:addChemRisksColumns (Toxicity Risks)');
      const compute: TemplateCompute = {descriptors: {enabled: true, args: ['MolWt', 'NumHDonors']},
        functions: [{name: toxicity.name, package: toxicity.package.name,
          args: {mutagenicity: false, tumorigenicity: true}, rerunOnOpen: true}], scripts: [], queries: []};
      let answered: IComputeDialogResult | undefined;
      const dialog = await chemFunctionsDialog(app, (r) => answered = r, () => {}, {compute}, false, true);
      await waitUntil(() => dialog.root.querySelector('.statistics-compute-dialog-descriptors-group') !== null,
        'the compute dialog never listed the descriptors', 30000);
      dialog.okProxy();
      const editor = computeEditorOnRoot({functions: app.computeFunctions, value: compute, descriptors: descriptorsTree,
        orderKey: HTFunctionOrderingLSKey, rerunOnOpen: rerunOnOpenOption});
      try {
        await editor.ready;
        expectSame(editor.read(), templateCompute(answered!), 'the editor\'s compute against the dialog\'s');
        expectSame(editor.read().descriptors, compute.descriptors, 'the descriptors');
      } finally {
        editor.dispose();
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Compute: a descriptor picked in the tree turns the descriptors on, and its group counts it', async () => {
    const sb = await Sandbox.create('editor');
    const host = document.body.appendChild(document.createElement('div'));
    host.style.cssText = 'position: fixed; left: 0; top: 0; width: 900px; height: 400px; z-index: 10000; background: white';
    try {
      const app = await sb.hitDesign();
      const editor = computeEditorOnRoot({functions: app.computeFunctions, descriptors: descriptorsTree,
        orderKey: HTFunctionOrderingLSKey, value: {descriptors: {enabled: false, args: []}, functions: []}});
      host.append(editor.root);
      try {
        await editor.ready;
        typeInto(editor.root.querySelector<HTMLInputElement>('[data-u2-name="descriptor-search"] input')!, 'MolWt');
        const molWt = () => Array.from(editor.root.querySelectorAll<HTMLElement>('.u2-tree-row'))
          .find((r) => r.querySelector('.u2-tree-label')?.textContent === 'MolWt');
        await waitUntil(() => molWt() !== undefined, 'the search never showed MolWt', 10000);
        molWt()!.querySelector<HTMLElement>('.u2-tree-check')!.click();
        await waitUntil(() => editor.read().descriptors.args.includes('MolWt'), 'MolWt was not picked', 5000);
        expect(editor.read().descriptors.enabled, true, 'the descriptors calculated');
        // the rows are flat: MolWt's group is the last group row above it, its shown descriptors the rows down to the next
        const label = (r: HTMLElement) => r.querySelector('.u2-tree-label')?.textContent ?? '';
        const isGroup = (r: HTMLElement) => / · \d+ of \d+$/.test(label(r));
        const rows = Array.from(editor.root.querySelectorAll<HTMLElement>('.u2-tree-row'));
        const groupAt = rows.slice(0, rows.indexOf(molWt()!)).map(isGroup).lastIndexOf(true);
        const next = rows.findIndex((r, i) => i > groupAt && isGroup(r));
        const shown = rows.slice(groupAt + 1, next < 0 ? rows.length : next).map(label);
        const groupName = label(rows[groupAt]).split(' · ')[0];
        const groupRow = () => Array.from(editor.root.querySelectorAll<HTMLElement>('.u2-tree-row'))
          .find((r) => label(r).startsWith(`${groupName} · `));
        await waitUntil(() => / · 1 of \d+$/.test(label(groupRow()!)), 'the group label did not count the pick', 5000);
        groupRow()!.querySelector<HTMLElement>('.u2-tree-check')!.click();
        await waitUntil(() => editor.read().descriptors.args.length === 0, 'the group\'s checkbox did not drop its pick', 5000);
        groupRow()!.querySelector<HTMLElement>('.u2-tree-check')!.click();
        await waitUntil(() => editor.read().descriptors.args.length > 0, 'the group\'s checkbox picked nothing', 5000);
        expectSame([...editor.read().descriptors.args].sort(), [...shown].sort(),
          'a group picks the descriptors the search shows, not all of its own');
      } finally {
        editor.dispose();
      }
    } finally {
      host.remove();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Compute: without descriptors to pick from there is no Descriptors row, and the descriptors stay as given', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const given = {enabled: true, args: ['MolWt']};
      const editor = computeEditorOnRoot({functions: app.computeFunctions, orderKey: HTFunctionOrderingLSKey,
        value: {descriptors: given, functions: []}});
      try {
        await editor.ready;
        const rows = Array.from(editor.root.querySelectorAll('[role="option"]')).map((r) => r.textContent);
        expect(rows.includes('Descriptors'), false, `a Descriptors row among ${rows.join(', ')}`);
        expect(rows.length > 0, true, 'the functions listed');
        expectSame(editor.read().descriptors, given, 'the descriptors');
      } finally {
        editor.dispose();
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Compute: the list is named "Functions", and a function\'s check says what it toggles', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitDesign();
      const editor = computeEditorOnRoot({functions: app.computeFunctions, orderKey: HTFunctionOrderingLSKey,
        value: noComputeTemplate('ZZB3 Named', 'ZZB3').compute});
      document.body.append(editor.root);
      try {
        expect(editor.root.querySelector('[role="listbox"]')?.getAttribute('aria-label'), 'Functions', 'the list\'s name');
        const box = editor.root.querySelector('[role="option"] input[type="checkbox"]');
        expect(box != null, true, 'a function\'s check');
        box!.dispatchEvent(new PointerEvent('pointerenter'));
        const tip = 'Toggle calculation of this function';
        await waitUntil(() => Array.from(document.querySelectorAll('[role="tooltip"]')).some((t) => t.textContent === tip),
          `the check's tooltip is not "${tip}"`, 5000);
        box!.dispatchEvent(new PointerEvent('pointerleave'));
      } finally {
        editor.root.remove();
        editor.dispose();
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Compute: a function hidden by "Reorder…" is not calculated, as in the compute dialog', async () => {
    const sb = await Sandbox.create('editor');
    const orderKey = `ZZB3 order ${DG.Utils.randomString(6)}`;
    try {
      const app = await sb.hitDesign();
      const f = app.computeFunctions.functions.find((f) => f.inputs.length >= 2 && f.inputs[0].propertyType === 'dataframe');
      if (!f)
        throw new Error('the stand has no compute function');
      localStorage.setItem(orderKey, JSON.stringify({order: [], hidden: [f.friendlyName ?? f.name]}));
      const compute: TemplateCompute = {descriptors: {enabled: false, args: []}, scripts: [], queries: [],
        functions: [{name: f.name, package: f.package.name, args: {}}]};
      const editor = computeEditorOnRoot({functions: app.computeFunctions, value: compute, descriptors: descriptorsTree,
        orderKey});
      try {
        await editor.ready;
        expectSame(editor.read().functions, [], 'the functions calculated');
      } finally {
        editor.dispose();
      }
    } finally {
      localStorage.removeItem(orderKey);
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Triage: CREATE is disabled while the name or the key is empty or taken', async () => {
    const sb = await Sandbox.create('editor');
    const taken = triageTemplate(`ZZB3 Taken ${DG.Utils.randomString(6)}`, `ZZBT${DG.Utils.randomString(6)}`);
    try {
      await _package.files.writeAsText(`${TRIAGE_TEMPLATES}/${taken.name}.json`, JSON.stringify(taken));
      const app = await sb.hitTriage();
      const editor = await createTemplateAccordeon(app, {});
      const create = createButton(editor.root);
      expect(disabled(create), true, 'CREATE with nothing filled in');
      typeTriage(editor.root, 'Name', taken.name);
      typeTriage(editor.root, 'Key', 'ZZBF');
      expect(disabled(create), true, 'CREATE with a taken name');
      typeTriage(editor.root, 'Name', 'ZZB3 Free');
      expect(disabled(create), false, 'CREATE with a free name and key');
      typeTriage(editor.root, 'Key', taken.key);
      expect(disabled(create), true, 'CREATE with a taken key');
    } finally {
      await removeTriageTemplate(taken.name);
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Triage, PepTriage: a clone proposes the key\'s next free number, as Hit Design does (KHL, KHL2 taken → KHL3)', async () => {
    const sb = await Sandbox.create('editor');
    const tag = DG.Utils.randomString(6);
    const key = `ZZK${tag}K`;
    const triage = [triageTemplate(`ZZB3 Clone ${tag}`, key), triageTemplate(`ZZB3 Clone ${tag} 2`, `${key}2`)];
    const keyOf = (root: HTMLElement) => root.querySelector<HTMLInputElement>('[name="input-host-Key"] input')!.value;
    try {
      for (const t of triage)
        await _package.files.writeAsText(`${TRIAGE_TEMPLATES}/${t.name}.json`, JSON.stringify(t));
      const ht = await createTemplateAccordeon(await sb.hitTriage(), {}, triage[0]);
      expect(keyOf(ht.root), `${key}3`, 'the Hit Triage clone\'s key');

      const peptides = [pepTriageTemplate(`ZZB3 Peptides ${tag}`, key), pepTriageTemplate(`ZZB3 Peptides ${tag} 2`, `${key}2`)];
      for (const t of peptides)
        await sb.writeTemplate(t);
      const pt = await createPepTriageTemplateAccordeon(await sb.pepTriage(), {}, peptides[0]);
      expect(keyOf(pt.root), `${key}3`, 'the PepTriage clone\'s key');
    } finally {
      for (const t of triage)
        await removeTriageTemplate(t.name);
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('Hit Triage, PepTriage: a clone proposes a name no template has, as Hit Design does ("(copy)" taken → "(copy 2)")', async () => {
    const sb = await Sandbox.create('editor');
    const tag = DG.Utils.randomString(6);
    const name = `ZZB3 Copy ${tag}`;
    const triage = [triageTemplate(name, `ZZC${tag}A`), triageTemplate(`${name} (copy)`, `ZZC${tag}B`)];
    const nameOf = (root: HTMLElement) => root.querySelector<HTMLInputElement>('[name="input-host-Name"] input')!.value;
    try {
      for (const t of triage)
        await _package.files.writeAsText(`${TRIAGE_TEMPLATES}/${t.name}.json`, JSON.stringify(t));
      const ht = await createTemplateAccordeon(await sb.hitTriage(), {}, triage[0]);
      expect(nameOf(ht.root), `${name} (copy 2)`, 'the Hit Triage clone\'s name');

      const peptides = [pepTriageTemplate(name, `ZZC${tag}A`), pepTriageTemplate(`${name} (copy)`, `ZZC${tag}B`)];
      for (const t of peptides)
        await sb.writeTemplate(t);
      const pt = await createPepTriageTemplateAccordeon(await sb.pepTriage(), {}, peptides[0]);
      expect(nameOf(pt.root), `${name} (copy 2)`, 'the PepTriage clone\'s name');
    } finally {
      for (const t of triage)
        await removeTriageTemplate(t.name);
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('Hit Triage: a template the share refuses is reported, not announced as created', async () => {
    const sb = await Sandbox.create('editor');
    try {
      const app = await sb.hitTriage();
      const editor = await createTemplateAccordeon(app, {});
      let created = false;
      void editor.template.then(() => created = true);
      typeTriage(editor.root, 'Name', unwritable('triage'));
      typeTriage(editor.root, 'Key', 'ZZBU');
      DG.Balloon.closeAll();
      createButton(editor.root).click();
      await waitUntil(() => balloonShown(saveFailed('triage')), 'no balloon said the save failed', 15000);
      expect(created, false, 'the editor reported a template created');
      expect(balloonShown('Template created successfully'), false, 'a success balloon');
    } finally {
      await removeTriageTemplate(unwritable('triage'));
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('PepTriage: CREATE waits for the name, the key and the sequence column', async () => {
    const sb = await Sandbox.create('editor');
    try {
      await sb.writeTemplate(pepTriageTemplate('ZZB3 Peptide Template', 'ZZBQ'));
      const app = await sb.pepTriage();
      const editor = await createPepTriageTemplateAccordeon(app, {});
      const create = createButton(editor.root);
      typeTriage(editor.root, 'Name', 'ZZB3 Peptides');
      typeTriage(editor.root, 'Key', 'ZZBP');
      expect(disabled(create), true, 'CREATE with no sequence column');
      typeTriage(editor.root, 'Sequence Column', 'Helm');
      expect(disabled(create), false, 'CREATE with a name, a key and a sequence column');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});
});
