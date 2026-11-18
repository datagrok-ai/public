/**
 * Category `Reactions: add panel`: the Reaction Sketcher's Add panel (utils/reactions/reaction-add-panel.ts): its
 * slots (Reactant, Reagent, Product, Next step) add what is typed by Enter or +, a name as its
 * agent's SMILES, and a text the drawing refuses is said and kept; the slots are Datagrok's inputs, names suggested in
 * the Reagent slot only, in Datagrok's list. On a stand-in drawing that records what it is given; the drawing's own
 * placing is Crux Sketch's (its addMolecule with a role); the app's Next step on the vendored sketcher.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {ADD_REFUSED_CLASS, ADD_SLOTS, AddPanel, AddRole, NAMES_NOT_KNOWN, REAGENT_HINT_NAMES, addedText, agentNames,
  nameNotKnown, placeList, readsAsName, suggestedNames} from '../utils/reactions/reaction-add-panel';
import {readsAs} from '../utils/reactions/reactions-app';
import {DOCUMENT_KEY, forgetDocument, lastProduct, sameReactionText} from '../utils/reactions/reaction-document';
import {loadCrux} from '../engine/crux-module';
import {loadCruxSketch} from '../crux/crux-sketcher';
import {PackageFunctions} from '../package';
import {shellPanelsKept} from './reactions-app-tests';

/** A small agent table in Crux's form: HATU over the arrow, DMF under it. */
const TABLE_HATU = 'CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F';
const TABLE = JSON.stringify({
  'CN(C)C=O': {name: 'DMF', place: 'below'},
  'CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F': {name: 'HATU', place: 'above'},
});

/** A stand-in drawing: what it was given, and the reason it refuses a text with `!` in it. */
function standIn(): {added: [string, AddRole][], add: (text: string, role: AddRole) => string | null} {
  const added: [string, AddRole][] = [];
  return {added, add: (text, role) => {
    if (text.includes('!'))
      return 'Not read as a molecule';
    added.push([text, role]);
    return null;
  }};
}

/** A key press in `input`. */
const press = (input: HTMLInputElement, key: string) =>
  input.dispatchEvent(new KeyboardEvent('keydown', {key, bubbles: true, cancelable: true}));

/** `text` typed into `input`, as a person's keys leave it (its value, then an input event). */
const typeInto = (input: HTMLInputElement, text: string) => {
  input.focus();
  input.value = text;
  input.dispatchEvent(new Event('input', {bubbles: true}));
};

/** The suggestions a slot's list shows (none when it is hidden or there is none). */
function shownSuggestions(panel: AddPanel, role: AddRole): string[] {
  const slot = panel.input(role).closest('.chem-reactions-add-slot')!;
  const list = slot.querySelector('.ui-input-list');
  if (!list || list.classList.contains('tt-hide') || getComputedStyle(list).display === 'none')
    return [];
  return Array.from(list.querySelectorAll('.tt-suggestion')).map((s) => s.textContent ?? '');
}

/** Reactant, Reagent and Product were raw monospace boxes sharing one native datalist: Chrome
 * drew its popup black, in the OS theme, and offered reagent names (K2CO3, HATU) on the Reactant and Product slots,
 * which take none. Now Datagrok's inputs: names in the Reagent slot only, in Datagrok's own list; a suggestion picked
 * fills the slot, Enter then adds it as its structure. */
async function datagrokInputs(): Promise<void> {
  const target = standIn();
  const panel = new AddPanel(target, () => agentNames(TABLE));
  document.body.append(panel.root);
  try {
    expect(panel.root.querySelector('datalist'), null, 'no native datalist: Chrome draws its popup in the OS theme');
    for (const s of ADD_SLOTS) {
      const input = panel.input(s.role);
      expect(input.getAttribute('list'), null, `${s.label}: no datalist`);
      expect(input.classList.contains('ui-input-editor'), true, `${s.label}: Datagrok's input (${input.className})`);
      expect(input.closest('.ui-form') !== null, true, `${s.label}: in a Datagrok form`);
      expect(/monospace/i.test(getComputedStyle(input).fontFamily), false,
        `${s.label}: Datagrok's font, not a monospace box (${getComputedStyle(input).fontFamily})`);
    }
    for (const role of ['reactant', 'product'] as AddRole[]) {
      typeInto(panel.input(role), 'HA');
      await delay(300);
      expect(shownSuggestions(panel, role).join(), '', `${role}: no names offered`);
      panel.input(role).value = '';
    }
    const agent = panel.input('agent');
    typeInto(agent, 'ha');
    await awaitCheck(() => shownSuggestions(panel, 'agent').join() === 'HATU',
      `the Reagent slot suggests HATU: [${shownSuggestions(panel, 'agent').join()}]`, 3000);
    const list = agent.closest('.chem-reactions-add-slot')!.querySelector('.ui-input-list') as HTMLElement;
    expect(getComputedStyle(list).backgroundColor !== 'rgb(0, 0, 0)', true, 'not black');
    // picked by the keys: the slot filled, nothing added until Enter
    press(agent, 'ArrowDown');
    press(agent, 'Enter');
    expect(agent.value, 'HATU', 'the suggestion picked fills the slot');
    expect(target.added.length, 0, 'picking adds nothing yet');
    press(agent, 'Enter');
    expect(JSON.stringify(target.added), JSON.stringify([[TABLE_HATU, 'agent']]), 'Enter adds its structure');
    expect(agent.value, '', 'emptied once added');
    // picked by the mouse
    typeInto(agent, 'dm');
    await awaitCheck(() => shownSuggestions(panel, 'agent').join() === 'DMF', 'DMF suggested', 3000);
    list.querySelector('.tt-suggestion')!.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true}));
    expect(agent.value, 'DMF', 'a click fills the slot');
    expect(suggestedNames('a', ['HATU', 'Pd/C', 'AcOH', 'DMAP']).join(), 'AcOH,HATU,DMAP',
      'names that start with the text first, then those it is in');
    // low in the window (the app's panel at 1440 x 900), the list opens above the slot, not under the status bar
    const slot = agent.closest('.chem-reactions-add-slot') as HTMLElement;
    Object.assign(panel.root.style, {position: 'fixed', left: '0', width: '600px', bottom: '4px'});
    placeList(slot, agent);
    expect(list.style.bottom, '100%', 'opens above, low in the window');
    Object.assign(panel.root.style, {bottom: '', top: '4px'});
    placeList(slot, agent);
    expect(list.style.bottom, '', 'opens below, high in the window');
    const scrollers = [panel.root, panel.root.querySelector('.ui-form') as HTMLElement];
    expect(scrollers.every((e) => e.scrollLeft === 0), true, 'nothing scrolled sideways by a slot');
  } finally {
    panel.root.remove();
  }
}

/**
 * The app's view, made by the bundle whose Crux Sketch defined <crux-sketch> on this page. The package tests run in
 * package-test.js beside package.js, each with its own copy of Crux Sketch, and the copy that defines the element first
 * is the only one that makes sketchers on the page: the other's createSketcher throws "Crux Sketch: <crux-sketch> is
 * defined by another script". Run alone, this category's Draw... test opens Chem's molecule sketcher (package.js)
 * before the app tests, so the test bundle's app never started its sketcher ("the sketcher did not start", in a run of
 * this category alone). A user has one bundle: the app opened cold starts its sketcher in 0.6 s, after a molecule
 * sketcher too (as measured in a browser). Importing this bundle's copy defines the element when no copy has; else the
 * app is opened through the platform, as a person opens it. The view is on the shell.
 */
async function openSketcherApp(): Promise<DG.View> {
  const ours = (await loadCruxSketch()).CruxSketchElement;
  const view = (customElements.get('crux-sketch') === ours ? await PackageFunctions.reactionsApp() :
    await grok.functions.call('Chem:reactionsApp')) as DG.View;
  if (!Array.from(grok.shell.views).some((v) => v.root === view.root))
    grok.shell.addView(view);
  return view;
}

/** The app's Next step adds the product of a new step after the route's last product, one Enter, the Text a route of
 * two steps (crux-sketch's addMolecule role 'next', API-075; red on crux-sketch caed19e8, which drops a role it does
 * not know and adds the molecule with none). */
async function nextStepAdded(): Promise<void> {
  forgetDocument();
  const view = await openSketcherApp();
  const text = () => (view.root.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement).value;
  const slot = () => view.root.querySelector('.chem-reactions-add-slot[data-role="next"]') as HTMLElement;
  const sketch = () => view.root.querySelector('.chem-reactions-app-sketcher crux-sketch') as any;
  try {
    await awaitCheck(() => sketch()?.document === 'reaction', 'the sketcher did not start', 15000);
    expect(slot() !== null && slot().style.display !== 'none', true, 'Next step offered');
    const box = view.root.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement;
    box.value = 'CCO>>CC=O';
    box.dispatchEvent(new Event('input', {bubbles: true}));
    box.dispatchEvent(new Event('blur'));
    await awaitCheck(() => sketch().isEmpty === false, 'the step not drawn', 5000);
    const input = slot().querySelector('input') as HTMLInputElement;
    input.value = 'CC(=O)O';
    input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true, cancelable: true}));
    expect(input.value, '', 'the slot emptied once added');
    const crux = await loadCrux('tests');
    const steps = () => {
      const scheme = new crux.ChemReactionScheme(text(), '');
      try {
        return scheme.get_steps();
      } finally {
        scheme.free();
      }
    };
    await awaitCheck(() => steps().length === 2, `the Text a route of two steps: "${text()}"`, 5000);
    const [one, two] = steps();
    expect(sameReactionText(crux, one, 'CCO>>CC=O') && sameReactionText(crux, two, 'CC=O>>CC(=O)O'), true,
      `the new step from the last product to the one typed: ${one} | ${two}`);
    expect(view.root.querySelector('.chem-reactions-app-status')?.textContent, readsAs(2), 'the status');
  } finally {
    view.close();
  }
}

/** A slot's Draw... button. */
const drawButton = (panel: AddPanel, role: AddRole) => Array.from(panel.input(role).closest('.chem-reactions-add-slot')!
  .querySelectorAll('button')).find((b) => b.textContent === 'Draw...') as HTMLButtonElement;

/** Whether the shell's current object is a molecule written as a MOL block (what Datagrok's top right label then
 * shows as raw text). */
const shellHoldsMolBlock = () => {
  const o = grok.shell.o;
  return o instanceof DG.SemanticValue && /M {2}END/.test(String(o.value ?? ''));
};

/** Next step's Draw... opened empty, the last product copied in by hand; the dialog's drawing (a paste)
 * became the shell's current object, its MOL block shown as raw text in the top right label across views. Now Next
 * step opens on the route's last product, the other slots empty, and the drawing stays the dialog's. */
async function drawOpensOnLastProduct(): Promise<void> {
  const target = {...standIn(), lastProduct: () => 'CC(=O)O'};
  const panel = new AddPanel(target, () => agentNames(null));
  document.body.append(panel.root);
  const crux = await loadCrux('tests');
  const before = ui.div([], 'chem-add-panel-test-current');
  grok.shell.o = before;
  try {
    drawButton(panel, 'next').click();
    await awaitCheck(() => panel.drawn !== null, 'Next step\'s Draw... opened no sketcher', 3000);
    const next = panel.drawn!;
    await next.sketcherReady();
    await awaitCheck(() => !next.isEmpty(), 'Next step\'s Draw... opened empty', 10000);
    expect(sameReactionText(crux, `${next.getSmiles()}>>C`, 'CC(=O)O>>C'), true,
      `it opens on the route's last product: ${next.getSmiles()}`);
    next.setSmiles('CC(=O)OC'); // as a paste changes it
    await delay(1000);
    expect(shellHoldsMolBlock(), false, 'the drawing is not the shell\'s current object');
    DG.Dialog.getOpenDialogs().find((d) => d.root.contains(next.root))?.close();
    await awaitCheck(() => panel.drawn === null, 'the dialog did not close', 3000);
    drawButton(panel, 'reactant').click();
    await awaitCheck(() => panel.drawn !== null, 'Reactant\'s Draw... opened no sketcher', 3000);
    const reactant = panel.drawn!;
    await reactant.sketcherReady();
    expect(reactant.isEmpty(), true, 'Reactant\'s Draw... opens empty');
    DG.Dialog.getOpenDialogs().find((d) => d.root.contains(reactant.root))?.close();
  } finally {
    const open = panel.drawn;
    DG.Dialog.getOpenDialogs().filter((d) => open && d.root.contains(open.root)).forEach((d) => d.close());
    panel.root.remove();
    if (grok.shell.o === before)
      grok.shell.o = null;
  }
}

/** The band at the canvas's foot that the readout ("Steps .. Reactants ..") covers, and the room kept inside the
 * canvas's edges (CSS px). */
const READOUT_BAND = 32;
const EDGE_ROOM = 4;

/** A step added to the default route by Next step landed under the readout at the canvas's foot, cut off. Every atom
 * of the drawing is in view after the addition: inside the canvas, clear of the readout's band (crux-sketch keeps what
 * addMolecule adds in view; red on crux-sketch fdb645b). */
async function nextStepInView(): Promise<void> {
  forgetDocument();
  const view = await openSketcherApp();
  const sketch = () => view.root.querySelector('.chem-reactions-app-sketcher crux-sketch') as any;
  const text = () => (view.root.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement).value;
  try {
    await awaitCheck(() => sketch()?.document === 'reaction' && sketch().isEmpty === false,
      'the sketcher did not show the default route', 15000);
    const before = text();
    const input = view.root.querySelector('.chem-reactions-add-slot[data-role="next"] input') as HTMLInputElement;
    input.value = 'CC(=O)OC';
    input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true, cancelable: true}));
    await awaitCheck(() => text() !== before, 'the step was not added', 5000);
    await awaitCheck(() => !sketch().isPending, 'the sketcher did not settle', 5000);
    const canvas = sketch().layout.canvas as {x: number, y: number, width: number, height: number};
    const atoms = (sketch().positions.atoms as ({x: number, y: number} | null)[]).filter((p) => p !== null);
    const out = atoms.filter((p) => p!.x < canvas.x + EDGE_ROOM || p!.x > canvas.x + canvas.width - EDGE_ROOM ||
      p!.y < canvas.y + EDGE_ROOM || p!.y > canvas.y + canvas.height - READOUT_BAND);
    expect(atoms.length > 0 && out.length === 0, true, `${out.length} of ${atoms.length} atoms out of view or under ` +
      `the readout: ${JSON.stringify(out.slice(0, 3))}, canvas ${JSON.stringify(canvas)}`);
  } finally {
    view.close();
  }
}

category('Reactions: add panel', () => {
  test('the agent table\'s names: those over the arrow first, any case read; a table that does not read gives none',
    async () => {
      const names = agentNames(TABLE);
      expect(names.names.join(), 'HATU,DMF');
      expect(addedText(' hatu ', names), 'CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F');
      expect(addedText('CC(=O)O', names), 'CC(=O)O', 'not a name: the text itself');
      expect(agentNames('{not json').names.length, 0);
      expect(agentNames(null).names.length, 0);
    });

  // an agent table entry may name the agent in other words too ("acetic anhydride" for Ac2O); the drawn name stays
  // the entry's
  test('the agent table\'s other names: offered after the names, any case read as the entry\'s SMILES; never in ' +
    'place of a name', async () => {
    const table = JSON.stringify({
      'CC(=O)OC(C)=O': {name: 'Ac2O', place: 'above', also: ['acetic anhydride', ' ', 7, 'DMF']},
      'CN(C)C=O': {name: 'DMF', place: 'below', also: 'dimethylformamide'},
    });
    const names = agentNames(table);
    expect(names.names.join(), 'Ac2O,DMF,acetic anhydride', 'the other names after the names; none blank or not text');
    expect(addedText('Acetic Anhydride', names), 'CC(=O)OC(C)=O', 'read as the entry\'s SMILES');
    expect(addedText('dmf', names), 'CN(C)C=O', 'an entry\'s name is its own, never another\'s other name');
    expect(suggestedNames('acet', names.names).join(), 'acetic anhydride', 'suggested as typed');
  });

  test('Chem\'s agent table names acetic anhydride as Ac2O (needs the vendored chem-wasm, ENGINE-7)', async () => {
    const crux = await loadCrux('tests') as unknown as {chemAgentTable?: () => string};
    const names = agentNames(crux.chemAgentTable?.() ?? null);
    expect(names.names.includes('Ac2O'), true, 'Ac2O in the table');
    expect(addedText('acetic anhydride', names), 'CC(=O)OC(C)=O',
      'acetic anhydride read as Ac2O\'s structure: needs the vendored chem-wasm (ENGINE-7)');
  });

  test('the slots: Enter or + adds what is typed as its role, a name as its structure; slots empty', async () => {
    const target = standIn();
    const panel = new AddPanel(target, () => agentNames(TABLE));
    document.body.append(panel.root);
    try {
      expect(ADD_SLOTS.map((s) => s.label).join(), 'Reactant,Reagent,Product,Next step');
      panel.input('reactant').value = 'CC(=O)O';
      press(panel.input('reactant'), 'Enter');
      panel.input('reactant').value = 'CCN';
      panel.addFrom('reactant'); // the + button's
      panel.input('agent').value = 'HATU';
      press(panel.input('agent'), 'Enter');
      panel.input('product').value = 'CCNC(C)=O';
      press(panel.input('product'), 'Enter');
      panel.input('next').value = 'CCN(C)C(C)=O';
      press(panel.input('next'), 'Enter');
      expect(JSON.stringify(target.added), JSON.stringify([['CC(=O)O', 'reactant'], ['CCN', 'reactant'],
        ['CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F', 'agent'], ['CCNC(C)=O', 'product'],
        ['CCN(C)C(C)=O', 'next']]));
      expect(ADD_SLOTS.every((s) => panel.input(s.role).value === ''), true, 'every slot emptied');
    } finally {
      panel.root.remove();
    }
  });

  test('a text the drawing refuses is said under the slots and kept in its slot; the next addition clears it',
    async () => {
      const target = standIn();
      const panel = new AddPanel(target, () => agentNames(null));
      document.body.append(panel.root);
      try {
        const line = panel.root.querySelector(`.${ADD_REFUSED_CLASS}`) as HTMLElement;
        panel.input('reactant').value = 'C!C';
        press(panel.input('reactant'), 'Enter');
        expect(line.style.display !== 'none' && line.textContent, 'Not read as a molecule');
        expect(panel.input('reactant').value, 'C!C', 'kept to be put right');
        panel.input('reactant').value = 'CC';
        press(panel.input('reactant'), 'Enter');
        expect(line.style.display, 'none', 'cleared');
        expect(target.added.length, 1);
      } finally {
        panel.root.remove();
      }
    });

  // "HATU", "DMF" and "acetic acid" were refused with 'unknown element symbol "H" at position 0', the
  // hint promising names; a refused text kept in its slot had the next entry appended to it ("HATUDMF")
  test('names: no hint promises them until known; a name refused is said as one; the next entry replaces it',
    async () => {
      const refusing = {add: (text: string) => /^[A-Za-z ]+$/.test(text) ? 'unknown element symbol "H" at position 0' :
        null};
      expect(readsAsName('HATU') && readsAsName('acetic acid') && readsAsName('Pd-C'), true, 'names');
      expect(readsAsName('C!C') || readsAsName('CC(=O)O') || readsAsName('[Pd]') || readsAsName('a\nb'), false,
        'not names');
      expect(ADD_SLOTS.some((s) => /name/i.test(s.hint)), false, 'no hint names names until they are known');
      for (const [table, said] of [[null, NAMES_NOT_KNOWN], [TABLE, nameNotKnown('Oxone')]] as const) {
        const panel = new AddPanel(refusing, () => agentNames(table));
        document.body.append(panel.root);
        try {
          const line = panel.root.querySelector(`.${ADD_REFUSED_CLASS}`) as HTMLElement;
          const slot = panel.input('agent');
          slot.dispatchEvent(new Event('focus'));
          expect(slot.placeholder === REAGENT_HINT_NAMES, table !== null, `the hint: ${slot.placeholder}`);
          slot.value = table === null ? 'HATU' : 'Oxone';
          press(slot, 'Enter');
          expect(line.style.display !== 'none' && line.textContent, said);
          expect(slot.value !== '' && slot.selectionStart === 0 && slot.selectionEnd === slot.value.length, true,
            'kept, selected: what is typed next replaces it');
          panel.input('reactant').value = 'acetic acid';
          press(panel.input('reactant'), 'Enter');
          expect(line.textContent, said.replace('Oxone', 'acetic acid'), 'a name in the Reactant slot');
        } finally {
          panel.root.remove();
        }
      }
    });

  test('Datagrok inputs: names offered in the Reagent slot only, in the platform\'s list, not a native datalist',
    datagrokInputs);

  test('Draw...: Next step opens on the route\'s last product, the others empty; the drawing is not the shell\'s ' +
    'current object', drawOpensOnLastProduct);

  test('the route\'s last product: the last step\'s last product, map numbers left out; none for no product',
    async () => {
      const crux = await loadCrux('tests');
      expect(lastProduct(crux, '[CH3:1][OH:2]>>[CH2:1]=[O:2]'), 'C=O');
      expect(lastProduct(crux, 'CCO>>CC=O.O--**--CC=O>>CC(=O)O'), 'CC(=O)O', 'a route: its last step\'s');
      expect(lastProduct(crux, 'CCO>>CC=O.[Na+].[Cl-]'), '[Cl-]', 'the last molecule of the products as written');
      expect(lastProduct(crux, 'CCO>>'), null, 'no product');
      expect(lastProduct(crux, ''), null, 'nothing');
    });

  test('the app\'s Next step adds a step after the last product, its product typed (crux-sketch from P12c)',
    async () => inTheApp(nextStepAdded));

  test('the app\'s Next step on the default route keeps the whole drawing in view, clear of the readout ' +
    '(crux-sketch from P12d)', async () => inTheApp(nextStepInView));
});

/** `run`, the app opened on a document of its own: the browser's document and the shell's panels put back after. */
async function inTheApp(run: () => Promise<void>): Promise<void> {
  const kept = window.localStorage.getItem(DOCUMENT_KEY);
  const putPanelsBack = shellPanelsKept();
  try {
    await run();
  } finally {
    if (kept === null)
      window.localStorage.removeItem(DOCUMENT_KEY);
    else
      window.localStorage.setItem(DOCUMENT_KEY, kept);
    putPanelsBack();
  }
}
