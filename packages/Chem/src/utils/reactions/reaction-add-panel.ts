/**
 * The Reaction Sketcher's Add panel: reactant, reagent and product slots taking drawn, pasted or typed molecules,
 * SMILES or a reagent's name, so that the app says how to add something. Reactant, Reagent and Product slots, and Next
 * step (the product of a new step after the route's last product); one line each:
 * - typed or pasted text (a SMILES, a MOL block, or, in the Reagent slot, an agent's name from Crux's agent table,
 *   which that slot alone suggests as it is typed) is added by Enter or the slot's + button;
 * - Draw... opens Chem's molecule sketcher in a dialog, and its OK adds what is drawn; for Next step it opens on the
 *   route's last product, to be changed into the next one (opened empty, the product had to be copied in by hand).
 *   The dialog's drawing is not made the shell's current object (a paste there put its raw MOL block in the
 *   top right label, kept across views).
 * The slots are Datagrok's own inputs in a form (not raw monospace boxes sharing a native datalist, whose popup Chrome
 * drew black in the OS theme, offering reagent names on the Reactant and Product slots too): the
 * Reagent slot is `ui.typeAhead` with a callback source, the names read as typed (the table may come after the panel),
 * its list Datagrok's (js-api src/widgets/inputs.ts TypeAhead; as UsageAnalysis's src/tabs/functions.ts uses it).
 * Each addition is one edit of the reaction drawing, which places it: a reactant after the reactants with a plus, a
 * reagent over the arrow, a product after the arrow, the arrow made first when there is none; a next step's product
 * past a new arrow after the last product, which becomes that step's reactant (Crux Sketch's addMolecule with a role:
 * one undo step). A name is added as its structure (the reaction text never holds a name or a
 * "*"; the publication look writes the name over the arrow). A text the drawing does not take is said under the slots,
 * and the slot keeps it to be put right.
 * What adds is behind {@link AddTarget}, so the panel is the same over any drawing that places a molecule by its role.
 */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import type {MoleculeRole} from '../../../vendor/crux-sketch/crux-sketch';

/** Where an added molecule goes in the reaction: Crux Sketch's addMolecule roles (`next`, the product of a new step
 * after the route's last product, crux-sketch API-075). */
export type AddRole = MoleculeRole;

/** What the panel adds to: `add` puts `text` into the reaction as `role`, one edit; null when it did, else why not.
 * `lastProduct`: the SMILES of the route's last product (the one a next step starts from), or null when there is
 * none. */
export interface AddTarget {
  add(text: string, role: AddRole): string | null;
  lastProduct?(): string | null;
}

/** The panel's class, and its slots' (`data-role` the slot's role): what the tests and bindings find them by. */
export const ADD_PANEL_CLASS = 'chem-reactions-add-panel';
export const ADD_SLOT_CLASS = 'chem-reactions-add-slot';
/** The class of the line under the slots that says why a text was not added. */
export const ADD_REFUSED_CLASS = 'chem-reactions-add-refused';
/** How many name suggestions the Reagent slot lists at most, and how long after a key it asks for them (ms). */
const SUGGESTIONS_SHOWN = 10;
const SUGGEST_MS = 50;

/** The slots, in the order shown: the role, the slot's label, the hint in its empty box, its tooltip. The Reagent
 * slot's hint names reagents only once Chem's Crux file has its agent table ({@link REAGENT_HINT_NAMES}; else the
 * hint would promise names the slot does not know). */
export const ADD_SLOTS: {role: AddRole, label: string, hint: string, tooltip: string}[] = [
  {role: 'reactant', label: 'Reactant', hint: 'SMILES, or paste a molecule',
    tooltip: 'Adds a reactant: after the reactants, with a plus'},
  {role: 'agent', label: 'Reagent', hint: 'SMILES, or paste a molecule',
    tooltip: 'Adds a reagent, catalyst or solvent over the arrow'},
  {role: 'product', label: 'Product', hint: 'SMILES, or paste a molecule',
    tooltip: 'Adds a product after the arrow'},
  // a route built step by step (crux-sketch's addMolecule role 'next', API-075)
  {role: 'next', label: 'Next step', hint: 'SMILES of the next product',
    tooltip: 'Adds a step: an arrow after the last product, then this molecule as the new step\'s product'},
];
/** The Reagent slot's hint once reagent names are known (the agent table). */
export const REAGENT_HINT_NAMES = 'SMILES, or a reagent\'s name (HATU, DMF)';
/** What the panel says of a name the drawing refused: none known here yet, or not this one. */
export const NAMES_NOT_KNOWN = 'Names are not read here yet: type the SMILES, or use Draw...';
export const nameNotKnown = (name: string) =>
  `"${name}" is not a reagent name Chem knows: type its SMILES, or use Draw...`;
/** The words of the panel's controls. */
const ADD_BUTTON_TOOLTIP = 'Add it (Enter)';
const DRAW_BUTTON = 'Draw...';
const DRAW_TOOLTIP = 'Draw the molecule, then OK adds it';

/** Agents by name: each name as written, and its SMILES under its lower-case name. */
export interface AgentNames {
  names: string[];
  smiles: Map<string, string>;
}

/** The agent table's names (Crux's `chemAgentTable()` JSON: `{"<SMILES>": {"name", "place", "also"?}}`), the ones
 * written over the arrow first; then the entries' other names (`also`: "acetic anhydride" for Ac2O),
 * each read as its entry's SMILES (the name drawn stays the entry's), never in place of an entry's name or of one
 * offered before it. None for a table that does not read. */
export function agentNames(tableJson: string | null): AgentNames {
  const out: AgentNames = {names: [], smiles: new Map()};
  let table: {[smiles: string]: {name?: string, place?: string, also?: unknown}} = {};
  try {
    table = tableJson ? JSON.parse(tableJson) : {};
  } catch {
    return out;
  }
  const entries = Object.entries(table).filter(([, e]) => typeof e?.name === 'string' && e.name);
  entries.sort(([, a], [, b]) => (a.place === 'below' ? 1 : 0) - (b.place === 'below' ? 1 : 0));
  const offer = (name: string, smiles: string) => {
    if (out.smiles.has(name.toLowerCase()))
      return;
    out.names.push(name);
    out.smiles.set(name.toLowerCase(), smiles);
  };
  for (const [smiles, e] of entries)
    offer(e.name!, smiles);
  for (const [smiles, e] of entries) {
    for (const other of Array.isArray(e.also) ? e.also : []) {
      if (typeof other === 'string' && other.trim())
        offer(other.trim(), smiles);
    }
  }
  return out;
}

/** Whether `typed` reads as a name (words on one line: letters, digits, spaces, hyphens, commas, apostrophes), not as
 * SMILES or a molecule's block. */
export function readsAsName(typed: string): boolean {
  return /^[A-Za-z][A-Za-z0-9 ,'-]*$/.test(typed.trim());
}

/** What `typed` adds: the SMILES of the agent it names (any case), else the text itself, trimmed. */
export function addedText(typed: string, names: AgentNames): string {
  const t = typed.trim();
  return names.smiles.get(t.toLowerCase()) ?? t;
}

/** The names that `typed` starts (any case) first, then those it is in, at most `limit`; none for an empty text. */
export function suggestedNames(typed: string, names: readonly string[], limit = SUGGESTIONS_SHOWN): string[] {
  const t = typed.trim().toLowerCase();
  if (!t)
    return [];
  const starts = names.filter((n) => n.toLowerCase().startsWith(t));
  const within = names.filter((n) => !n.toLowerCase().startsWith(t) && n.toLowerCase().includes(t));
  return [...starts, ...within].slice(0, limit);
}

/** The room a suggestion list takes (CSS px): {@link SUGGESTIONS_SHOWN} rows of Datagrok's list. */
const LIST_ROOM = 260;

/** The slot `root`'s suggestion list placed under `input`, or above it when the window has no room below and more above
 * (Datagrok's list hangs under its input). */
export function placeList(root: HTMLElement, input: HTMLInputElement): void {
  const list = root.querySelector('.ui-input-list') as HTMLElement | null;
  if (!list)
    return;
  const r = input.getBoundingClientRect();
  const above = window.innerHeight - r.bottom < LIST_ROOM && r.top > window.innerHeight - r.bottom;
  Object.assign(list.style, above ? {top: 'auto', bottom: '100%'} : {top: '', bottom: ''});
}

/** The Add panel over `target`, its names from `names` (read when a slot is first used: the table may come later). */
export class AddPanel {
  readonly root: HTMLElement;
  private readonly refused = ui.divText('', ADD_REFUSED_CLASS);
  private readonly inputs = new Map<AddRole, HTMLInputElement>();
  private known: AgentNames | null = null;
  /** The Draw dialog's sketcher while the dialog is open (tests). */
  drawn: DG.chem.Sketcher | null = null;

  constructor(private readonly target: AddTarget, private readonly names: () => AgentNames) {
    const fields = ADD_SLOTS.map((slot) => {
      // names in the Reagent slot only: the slot that takes an agent's name
      const field: DG.InputBase = slot.role === 'agent' ?
        ui.typeAhead(slot.label, {source: async (typed: string) => suggestedNames(typed, this.agentNamesNow().names),
          minLength: 1, limit: SUGGESTIONS_SHOWN, debounceRemote: SUGGEST_MS}) :
        ui.input.string(slot.label, {value: ''});
      const input = field.input as HTMLInputElement;
      Object.assign(input, {placeholder: slot.hint, spellcheck: false, autocomplete: 'off'});
      input.setAttribute('aria-label', slot.label);
      if (slot.role === 'agent') {
        input.addEventListener('focus', () => {
          if (this.agentNamesNow().names.length > 0) // the agent table is in Chem's Crux file
            input.placeholder = REAGENT_HINT_NAMES;
        });
        // the list opens above the slot when the window has no room for it below (the panel sits low in the app,
        // where the list went under the window's status bar at 1440 x 900)
        input.addEventListener('input', () => placeList(field.root, input));
      }
      input.addEventListener('keydown', (e) => {
        // a suggestion picked by Enter fills the slot (the list's own Enter, default prevented); the next Enter adds
        if (e.key === 'Enter' && !e.defaultPrevented) {
          e.preventDefault();
          this.addFrom(slot.role);
        }
      });
      this.inputs.set(slot.role, input);
      field.setTooltip(slot.tooltip);
      const add = ui.iconFA('plus', () => this.addFrom(slot.role), ADD_BUTTON_TOOLTIP);
      const draw = ui.button(DRAW_BUTTON, () => this.draw(slot.role), DRAW_TOOLTIP);
      field.addOptions(ui.divH([add, draw], {style: {alignItems: 'center', gap: '6px'}}));
      field.root.classList.add(ADD_SLOT_CLASS);
      field.root.dataset.role = slot.role;
      return field;
    });
    Object.assign(this.refused.style, {display: 'none', color: 'var(--red-3)', fontSize: '12px'});
    // the form's own overflow-x is hidden, a box a focused slot scrolls sideways (the labels cut):
    // clipped here instead, which nothing scrolls, the list free to show below or above
    const form = ui.form(fields);
    Object.assign(form.style, {overflowX: 'visible', overflowY: 'visible'});
    this.root = ui.divV([form, this.refused], ADD_PANEL_CLASS);
    Object.assign(this.root.style, {margin: '6px 0', overflowX: 'clip', overflowY: 'visible'});
  }

  /** The slot of `role`'s box (tests). */
  input(role: AddRole): HTMLInputElement {
    return this.inputs.get(role)!;
  }

  private agentNamesNow(): AgentNames {
    if (this.known === null || this.known.names.length === 0)
      this.known = this.names();
    return this.known;
  }

  /** The slot of `role`'s text added: emptied when the drawing took it, kept with the reason said when not. */
  addFrom(role: AddRole): void {
    const input = this.inputs.get(role)!;
    if (!input.value.trim())
      return;
    const names = this.agentNamesNow();
    const text = addedText(input.value, names);
    // a name refused: said as a name not known, not with the drawing's reason for a structure
    const asName = text === input.value.trim() && readsAsName(text) ?
      (names.names.length ? nameNotKnown(text) : NAMES_NOT_KNOWN) : null;
    if (this.add(text, role, asName))
      input.value = '';
    else
      input.select(); // kept to be put right, and what is typed next replaces it (not "HATUDMF")
  }

  // `text` added as `role`; whether it was (else the reason is said under the slots: `asName` in place of the
  // drawing's when given).
  private add(text: string, role: AddRole, asName: string | null = null): boolean {
    let why: string | null;
    try {
      why = this.target.add(text, role);
    } catch (e) {
      console.error('Chem | the Add panel could not add', e);
      why = 'That could not be added';
    }
    if (why !== null && asName !== null)
      why = asName;
    this.refused.textContent = why ?? '';
    this.refused.style.display = why ? '' : 'none';
    return why === null;
  }

  // Chem's molecule sketcher in a dialog, for Next step on the route's last product; its OK adds what is drawn as
  // `role`.
  private draw(role: AddRole): void {
    const sketcher = new DG.chem.Sketcher();
    // the drawing stays the dialog's: not the shell's current object, whose label showed a pasted MOL block as text
    sketcher.syncCurrentObject = false;
    const last = role === 'next' ? this.lastProduct() : null;
    if (last)
      sketcher.setSmiles(last);
    const label = ADD_SLOTS.find((s) => s.role === role)!.label;
    const dialog = ui.dialog(role === 'next' ? 'Draw the next step\'s product' : `Draw a ${label.toLowerCase()}`)
      .add(sketcher)
      .onOK(() => {
        const molfile = sketcher.getMolFile();
        if (molfile && !sketcher.isEmpty())
          this.add(molfile, role);
      })
      .show();
    this.drawn = sketcher;
    dialog.onClose.subscribe(() => {
      if (this.drawn === sketcher)
        this.drawn = null;
    });
  }

  // The route's last product, or null when there is none or it cannot be said.
  private lastProduct(): string | null {
    try {
      return this.target.lastProduct?.() ?? null;
    } catch (e) {
      console.warn('Chem | the route\'s last product was not read', e);
      return null;
    }
  }
}
