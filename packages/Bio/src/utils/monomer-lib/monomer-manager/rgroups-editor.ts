/* The R groups of the monomer in the editor: a chip per R label the structure carries — its cap, red while
   it cannot be saved. A click opens the chip's menu: the common caps, each with what it is and the SMILES
   Save writes, the cap its definition names when that is none of them, and Edit definition, which opens
   the fields a library stores under the chips. */
import {RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {Control, Form, TextInput, Menu, Scope, signal, computed, each, span, div, iconButton, icon, on, effect,
  Tooltip, Signal, ReadonlySignal, batch} from '@datagrok-libraries/u2';

import {CapGroups, COMMON_CAP_GROUPS, CapInference, Canonicalize} from './cap-groups';
import {STANDRARD_R_GROUPS} from './const';

export interface RGroupsEditorOptions {
  canonicalize?: Canonicalize;
  /** Whether a cap SMILES is valid; without it, any non-empty one is. */
  isValidSmiles?: (smiles: string) => boolean;
  /** The user changed an R group — a pick or one of the definition's fields. */
  onInput?: () => void;
}

type Field = 'smiles' | 'name' | 'alternateId';

const FIELD_KEYS: Record<Field, string> = {
  smiles: HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE,
  name: HELM_RGROUP_FIELDS.CAP_GROUP_NAME,
  alternateId: HELM_RGROUP_FIELDS.ALTERNATE_ID,
};

const byLabel = (a: string, b: string) => (CapGroups.numberOf(a) ?? 0) - (CapGroups.numberOf(b) ?? 0);

export class RGroupsEditor extends Control {
  /** The definitions in label order, as Save writes them: a common cap in its canonical form, any
   * other definition as it stands. */
  readonly rgroups: ReadonlySignal<RGroup[]>;
  /** The first R group that cannot be saved and why, or null. */
  readonly problem: ReadonlySignal<string | null>;

  /** The R group whose definition shows under the chips, or null. */
  readonly open = signal<string | null>(null);

  private readonly _labels = signal<string[]>([]);
  private readonly _defs = new Map<string, Signal<RGroup>>();
  /** Bumped by every monomer loaded: its chips are built afresh, while a structure edit keeps them. */
  private _generation = 0;
  /** What is wrong with the structure's R groups themselves (an R label twice), over the chips. */
  private readonly _structure = span('', 'monomer-rgroups-structure');

  constructor(private readonly _options: RGroupsEditorOptions = {}) {
    super();
    this.root.classList.add('monomer-rgroups-editor');
    this.rgroups = computed(() => this._labels.value.map((label) => this._saved(label, this._defs.get(label)!.value)));
    this.problem = computed(() => {
      for (const label of this._labels.value) {
        const problem = this._problemOf(label, this._defs.get(label)!.value);
        if (problem !== null)
          return problem;
      }
      return null;
    });
    this._structure.dataset.u2Part = 'error';
    this._structure.hidden = true;
    const chips = div([], 'monomer-rgroups-chips');
    const definitions = div([], 'monomer-rgroups-definitions');
    this.root.append(this._structure, chips, definitions);
    const empty = () =>
      span('Draw the monomer with its R groups (R1, R2, ...) to define their caps', 'monomer-rgroups-empty');
    const key = (label: string) => `${this._generation}:${label}`;
    this.runInScope(() => {
      each(chips, this._labels, key, (label) => this._chip(label), {empty});
      each(definitions, this._labels, key, (label) => this._definition(label));
    });
    // an R group the structure lost takes its definition with it
    this.effect(() => {
      const open = this.open.value;
      if (open !== null && !this._labels.value.includes(open))
        this.open.value = null;
    });
  }

  /** Shows the R groups of a monomer as its library defines them. */
  setDefinitions(rgroups: RGroup[]): void {
    batch(() => {
      this._generation++;
      this.open.value = null;
      this._defs.clear();
      for (const rg of rgroups) {
        const label = CapGroups.labelOf(rg);
        if (label !== '')
          this._defs.set(label, signal(rg));
      }
      this._labels.value = [...this._defs.keys()].sort(byLabel);
    });
  }

  /** Follows the R labels of the drawn structure: kept labels keep their definitions, a new one
   * starts as the usual cap for its position (R1 H, R2 OH, then H). */
  setLabels(labels: string[]): void {
    batch(() => {
      for (const label of labels) {
        if (!this._defs.has(label))
          this._defs.set(label, signal(RGroupsEditor.defaultFor(label)));
      }
      this._labels.value = [...labels].sort(byLabel);
    });
  }

  /** Says, over the chips, what is wrong with the R groups of the drawn structure; null clears it. */
  setStructureProblem(message: string | null): void {
    this._structure.textContent = message ?? '';
    this._structure.hidden = message === null;
  }

  get labels(): string[] {
    return this._labels.peek();
  }

  /** Sets R`n`'s cap by its name (`OH`) or by its four fields — the way code and copilots edit it. */
  setCap(label: string, cap: string | RGroup): void {
    const def = this._defs.get(label);
    if (def === undefined)
      throw new Error(`The monomer has no ${label}; its R groups are ${this._labels.peek().join(', ') || 'none'}`);
    const n = CapGroups.numberOf(label) ?? 1;
    const known = typeof cap === 'string' ? CapGroups.byName(cap) : null;
    if (typeof cap === 'string' && known === null)
      throw new Error(`"${cap}" is not one of the common caps: ${COMMON_CAP_GROUPS.map((c) => c.name).join(', ')}`);
    def.value = known !== null ? CapGroups.rgroup(known, n) :
      typeof cap === 'string' ? CapGroups.undefinedRGroup(label) : {...cap, [HELM_RGROUP_FIELDS.LABEL]: label};
  }

  static defaultFor(label: string): RGroup {
    const n = CapGroups.numberOf(label) ?? 1;
    const standard = STANDRARD_R_GROUPS.find((rg) => rg.label === label);
    const cap = (standard ? CapGroups.infer(standard)?.cap : null) ?? COMMON_CAP_GROUPS[0];
    return CapGroups.rgroup(cap, n);
  }

  private _infer(def: RGroup): CapInference | null {
    return CapGroups.infer(def, this._options.canonicalize);
  }

  private _saved(label: string, def: RGroup): RGroup {
    const inferred = this._infer(def);
    return inferred?.known ? CapGroups.rgroup(inferred.cap, CapGroups.numberOf(label) ?? 1) :
      {...def, [HELM_RGROUP_FIELDS.LABEL]: label};
  }

  private _problemOf(label: string, def: RGroup): string | null {
    const inferred = this._infer(def);
    if (inferred === null)
      return `${label}: pick its cap group, or edit its definition`;
    if (inferred.known)
      return null;
    if (CapGroups.isPlaceholder(inferred.cap))
      return `${label}: X is a placeholder, not a cap group — pick the cap ${label} gets (usually H)`;
    const smiles = CapGroups.smilesOf(def);
    if (smiles === '')
      return `${label}: the cap group ${inferred.cap.name} needs its SMILES`;
    if (this._options.isValidSmiles && !this._options.isValidSmiles(smiles))
      return `${label}: the cap SMILES ${smiles} is not valid`;
    return null;
  }

  /** What the chip and its menu say of a definition: the cap, what it is (its friendly name) and the SMILES
   * Save writes. */
  private _describe(label: string, def: RGroup): {cap: string, what: string, smiles: string} {
    const n = CapGroups.numberOf(label) ?? 1;
    const inference = this._infer(def);
    const smiles = CapGroups.smilesOf(def);
    if (inference === null) {
      const raw = smiles || CapGroups.nameOf(def);
      return {cap: raw ? '?' : '–', what: raw ? `not recognized: ${raw}` : 'not defined', smiles: ''};
    }
    return inference.known ?
      {cap: inference.cap.name, what: inference.cap.description + (inference.note ? ` (${inference.note})` : ''),
        smiles: CapGroups.numbered(inference.cap.smiles, n)} :
      {cap: inference.cap.name, what: inference.note ?? 'as defined', smiles: smiles || 'no SMILES'};
  }

  /** Why R`n` cannot be saved, without the label the chip already shows. */
  private _problemHere(label: string, def: RGroup): string | null {
    return this._problemOf(label, def)?.replace(`${label}: `, '') ?? null;
  }

  /** R1 and its cap, small: a click opens its menu, a second one closes it. */
  private _chip(label: string): HTMLElement {
    const def = this._defs.get(label)!;
    // the definition the menu brings back: the latest one that named a cap outside the common ones
    const own = signal<RGroup | null>(null);
    effect(() => {
      const inference = this._infer(def.value);
      if (inference !== null && !inference.known)
        own.value = def.peek();
    });
    const cap = span('', 'monomer-rgroup-chip-cap');
    const chip = document.createElement('button');
    chip.type = 'button';
    chip.className = 'monomer-rgroup-chip';
    chip.dataset.u2Item = label;
    chip.setAttribute('aria-haspopup', 'menu');
    chip.setAttribute('aria-expanded', 'false');
    const caret = icon('chevron-down', {cls: 'monomer-rgroup-chip-caret'});
    caret.setAttribute('aria-hidden', 'true');
    // the space reads "R4 X" to a screen reader; a flex row does not render it
    chip.append(span(label, 'monomer-rgroup-chip-label'), ' ', cap, caret);
    const said = computed(() => ({...this._describe(label, def.value), problem: this._problemHere(label, def.value)}));
    effect(() => {
      const {cap: name, what, smiles, problem} = said.value;
      cap.textContent = name;
      chip.classList.toggle('monomer-rgroup-chip-invalid', problem !== null);
      chip.setAttribute('aria-invalid', String(problem !== null));
      chip.classList.toggle('monomer-rgroup-chip-editing', this.open.value === label);
      // what the tooltip shows the pointer, said to the keyboard and a screen reader too
      chip.setAttribute('aria-description', problem ?? (smiles ? `${what} · ${smiles}` : what));
    });
    Tooltip.bind(chip, () => `${label}: ${chip.getAttribute('aria-description') ?? ''}`);

    let menu: Menu | null = null;
    const toggle = () => {
      if (menu !== null) {
        menu.close();
        return;
      }
      const shown = this._menu(label, own.peek());
      menu = shown;
      shown.show({anchor: chip});
      // aria-expanded follows the menu, and the watch goes with it, as a dropDownButton's does
      const watch = new Scope(undefined, {label: 'rgroup-menu'});
      watch.effect(() => {
        const open = shown.isOpen.value;
        chip.setAttribute('aria-expanded', String(open));
        if (open)
          return;
        menu = null;
        watch.dispose();
      });
    };
    on(chip, 'click', toggle);
    on(chip, 'keydown', (e) => {
      if ((e as KeyboardEvent).key !== 'ArrowDown' || menu !== null)
        return;
      e.preventDefault();
      toggle();
    });
    return chip;
  }

  /** The chip's menu, built on every open: the common caps with what each is and its SMILES, the current one
   * checked; the cap the definition names when it is none of them; and the definition's fields. */
  private _menu(label: string, own: RGroup | null): Menu {
    const def = this._defs.get(label)!;
    const n = CapGroups.numberOf(label) ?? 1;
    const current = this._infer(def.peek());
    const pick = (rg: RGroup) => {
      def.value = rg;
      this._options.onInput?.();
    };
    const menu = Menu.popup();
    for (const cap of COMMON_CAP_GROUPS) {
      menu.item(cap.name, () => pick(CapGroups.rgroup(cap, n)), {check: current?.known === true &&
        current.cap.name === cap.name, detail: `${cap.description} · ${CapGroups.numbered(cap.smiles, n)}`});
    }
    const ownCap = own === null ? null : this._infer(own)?.cap ?? null;
    if (own !== null && ownCap !== null) {
      const {what, smiles} = this._describe(label, own);
      // a cap of its own that borrows a common name (H over other chemistry) must not read as that one
      menu.item(CapGroups.byName(ownCap.name) === null ? ownCap.name : `${ownCap.name} (custom)`, () => pick(own),
        {check: current !== null && !current.known, detail: `${what} · ${smiles}`});
    }
    menu.separator();
    const editing = this.open.peek() === label;
    menu.item(editing ? 'Hide definition' : 'Edit definition…', () => this.open.value = editing ? null : label);
    return menu;
  }

  /** The fields a library stores for R`n`, under the chips while its definition is open; they show what
   * Save writes — a cap read from the alternate ID alone shows its SMILES and name too — and an edit
   * starts from that. */
  private _definition(label: string): HTMLElement {
    const def = this._defs.get(label)!;
    const saved = computed(() => this._saved(label, def.value));
    const fields = new Form({layout: 'normal'});
    const fieldInput = (caption: string, field: Field, read: (rg: RGroup) => string) => {
      const input = new TextInput({label: caption, name: `${label} ${caption}`, required: true, commitOn: 'change',
        onInput: (value) => {
          def.value = {...saved.peek(), [FIELD_KEYS[field]]: value.trim()};
          this._options.onInput?.();
        }});
      input.follow(() => read(saved.value));
      return input;
    };
    fields.addAll([
      fieldInput('Cap SMILES', 'smiles', CapGroups.smilesOf),
      fieldInput('Cap name', 'name', CapGroups.nameOf),
      fieldInput('Alternate ID', 'alternateId', CapGroups.alternateIdOf),
    ]);
    const close = iconButton('times', () => this.open.value = null, {tooltip: `Close the ${label} definition`});
    const head = div([span(`${label} definition`, 'monomer-rgroup-definition-title'), close],
      'monomer-rgroup-definition-head');
    const box = div([head, fields.root], 'monomer-rgroup-definition');
    box.id = `${this.rootId}-${label}`;
    box.dataset.u2Item = label;
    fields.effect(() => box.hidden = this.open.value !== label);
    return box;
  }
}
