/* The R groups of the monomer in the editor: a chip per R label the structure carries — its cap, red while
   it cannot be saved — and, above the chips, the row of the one opened: a choice of the common caps (and of
   the cap its definition names when that is none of them), a line saying what the chosen cap is, and a
   pencil that opens the four fields a library stores. */
import {RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {Control, Form, ChoiceInput, TextInput, ChoiceItem, LabelColumn, signal, computed, each, span, div,
  iconButton, icon, on, effect, Tooltip, Signal, ReadonlySignal, batch} from '@datagrok-libraries/u2';

import {CapGroup, CapGroups, COMMON_CAP_GROUPS, CapInference, Canonicalize} from './cap-groups';
import {STANDRARD_R_GROUPS} from './const';

/** The choice that stands for the cap the definition names when it is none of the common ones. */
const OWN = '~own';

export interface RGroupsEditorOptions {
  canonicalize?: Canonicalize;
  /** Whether a cap SMILES is valid; without it, any non-empty one is. */
  isValidSmiles?: (smiles: string) => boolean;
  /** The label column of the form around the rows; without one, R1, R2 get a narrow column of their own. */
  labelColumn?: LabelColumn;
  /** The user changed an R group — a pick or one of the four fields. */
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

  /** The R group whose row shows above the chips, or null. */
  readonly open = signal<string | null>(null);

  private readonly _labels = signal<string[]>([]);
  private readonly _defs = new Map<string, Signal<RGroup>>();
  /** Bumped by every monomer loaded: its rows are built afresh, while a structure edit keeps them. */
  private _generation = 0;
  private readonly _form: Form;
  /** What is wrong with the structure's R groups themselves (an R label twice), over the rows. */
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
    // R1, R2 alone want a narrow column, not the 100 px a form gives by default
    this._form = this.runInScope(() => new Form(_options.labelColumn !== undefined ?
      {layout: 'normal', labelColumn: _options.labelColumn} : {layout: 'normal', labelWidth: 28}));
    const rows = div([], 'monomer-rgroups-rows');
    this._form.addElement(rows);
    this._structure.dataset.u2Part = 'error';
    this._structure.hidden = true;
    const chips = div([], 'monomer-rgroups-chips');
    // the open row sits above the chips: with the fields growing upward, the chips stay where they were clicked
    this.root.append(this._structure, this._form.root, chips);
    const empty = () =>
      span('Draw the monomer with its R groups (R1, R2, ...) to define their caps', 'monomer-rgroups-empty');
    const key = (label: string) => `${this._generation}:${label}`;
    this.runInScope(() => {
      each(chips, this._labels, key, (label) => this._chip(label), {empty});
      each(rows, this._labels, key, (label) => this._row(label));
    });
    // an R group the structure lost takes its row with it
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

  /** Says, over the rows, what is wrong with the R groups of the drawn structure; null clears it. */
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
      return `${label}: pick its cap group, or define it with the pencil`;
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

  /** What the chip and the row say of a definition: the cap, what it is (its friendly name) and the SMILES
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

  /** Why R`n` cannot be saved, without the label the row and the chip already show. */
  private _problemHere(label: string, def: RGroup): string | null {
    return this._problemOf(label, def)?.replace(`${label}: `, '') ?? null;
  }

  /** R1 and its cap, small: opens the R group's row above the chips, or closes it. */
  private _chip(label: string): HTMLElement {
    const def = this._defs.get(label)!;
    const cap = span('', 'monomer-rgroup-chip-cap');
    const chip = document.createElement('button');
    chip.type = 'button';
    chip.className = 'monomer-rgroup-chip';
    chip.dataset.u2Item = label;
    chip.setAttribute('aria-controls', `${this.rootId}-${label}`);
    const caret = icon('chevron-down', {cls: 'monomer-rgroup-chip-caret'});
    caret.setAttribute('aria-hidden', 'true');
    // the space reads "R4 X" to a screen reader; a flex row does not render it
    chip.append(span(label, 'monomer-rgroup-chip-label'), ' ', cap, caret);
    const said = computed(() => ({...this._describe(label, def.value), problem: this._problemHere(label, def.value)}));
    const about = () => {
      const {what, smiles, problem} = said.value;
      return problem ?? (smiles ? `${what} · ${smiles}` : what);
    };
    effect(() => {
      const {cap: name, problem} = said.value;
      cap.textContent = name;
      chip.classList.toggle('monomer-rgroup-chip-invalid', problem !== null);
      chip.setAttribute('aria-expanded', String(this.open.value === label));
      // what the tooltip shows the pointer, said to the keyboard and a screen reader too
      chip.setAttribute('aria-description', about());
    });
    Tooltip.bind(chip, () => `${label}: ${chip.getAttribute('aria-description') ?? ''}`);
    on(chip, 'click', () => this.open.value = this.open.peek() === label ? null : label);
    return chip;
  }

  private _row(label: string): HTMLElement {
    const def = this._defs.get(label)!;
    const n = CapGroups.numberOf(label) ?? 1;
    const inferred = computed(() => this._infer(def.value));
    // the definition the "own" choice brings back: the latest one that named a cap outside the list
    const own = signal<RGroup | null>(null);
    const ownCap = computed(() => own.value === null ? null : this._infer(own.value)?.cap ?? null);
    // a cap of its own that borrows a common name (H over other chemistry) must not read as that one
    const ownLabel = (cap: CapGroup) => CapGroups.byName(cap.name) === null ? cap.name : `${cap.name} (custom)`;
    const items = computed<ChoiceItem[]>(() => [
      ...COMMON_CAP_GROUPS.map((c) => ({value: c.name, label: c.name})),
      ...ownCap.value === null ? [] : [{value: OWN, label: ownLabel(ownCap.value)}],
    ]);
    const editing = signal(false);

    const choice = new ChoiceInput({label, name: label, items, nullable: true, width: 110,
      tooltipText: `The cap group that closes ${label} when the monomer ends a chain`,
      onInput: (value) => {
        if (value === OWN && own.peek() !== null)
          def.value = own.peek()!;
        else if (value === null)
          def.value = CapGroups.undefinedRGroup(label);
        else {
          const cap = CapGroups.byName(value);
          if (cap !== null)
            def.value = CapGroups.rgroup(cap, n);
        }
        this._options.onInput?.();
      }});
    // every effect of the row is its controls': they go with the row, not with the editor
    choice.effect(() => {
      const inference = inferred.value;
      if (inference !== null && !inference.known)
        own.value = def.peek();
    });
    choice.follow(() => {
      const inference = inferred.value;
      return inference === null ? null : inference.known ? inference.cap.name : OWN;
    });
    choice.addValidator(() => {
      const problem = this._problemHere(label, def.peek());
      return problem === null ? null : problem.charAt(0).toUpperCase() + problem.slice(1);
    });
    choice.effect(() => {
      def.value;
      choice.revalidate();
    });

    // the friendly name, then the SMILES, which a narrow dock drops first
    const what = span('', 'monomer-rgroup-what');
    const smiles = span('', 'monomer-rgroup-smiles');
    const info = div([what, smiles], 'monomer-rgroup-info');
    info.dataset.u2Part = 'description';
    choice.effect(() => {
      const described = this._describe(label, def.value);
      what.textContent = described.what;
      smiles.textContent = described.smiles ? ` · ${described.smiles}` : '';
      info.title = info.textContent ?? '';
    });
    const pencil = iconButton('pencil', () => {}, {tooltip: `Edit the ${label} definition`, toggle: editing});
    choice.addOptions(info);
    choice.addOptions(pencil);

    const fields = new Form({layout: 'normal'});
    fields.root.classList.add('monomer-rgroup-fields');
    // the fields show what Save writes — a cap read from the alternate ID alone shows its SMILES and name
    // too — and an edit starts from that
    const saved = computed(() => this._saved(label, def.value));
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
      new TextInput({label: 'Label', name: `${label} Label`, value: label, readOnly: true}),
    ]);
    fields.effect(() => fields.root.hidden = !editing.value);

    const row = div([], 'monomer-rgroup-row');
    row.id = `${this.rootId}-${label}`;
    row.dataset.u2Item = label;
    // the pencil's fields open above its row, which stays where the pencil was clicked
    row.append(fields.root);
    this._form.add(choice, row);
    choice.effect(() => row.hidden = this.open.value !== label);
    return row;
  }
}
