/* The R groups of the monomer in the editor: a row per R label the structure carries, each a choice
   of the common caps (and of the cap its definition names when that is none of them), a line saying
   what the chosen cap is, and a pencil that opens the four fields a library stores. */
import {RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {Control, Form, ChoiceInput, TextInput, ChoiceItem, signal, computed, each, span, div, iconButton,
  Signal, ReadonlySignal, batch} from '@datagrok-libraries/u2';

import {CapGroup, CapGroups, COMMON_CAP_GROUPS, CapInference, Canonicalize} from './cap-groups';
import {STANDRARD_R_GROUPS} from './const';

/** The choice that stands for the cap the definition names when it is none of the common ones. */
const OWN = '~own';

export interface RGroupsEditorOptions {
  canonicalize?: Canonicalize;
  /** Whether a cap SMILES is valid; without it, any non-empty one is. */
  isValidSmiles?: (smiles: string) => boolean;
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
    this._form = this.runInScope(() => new Form({layout: 'normal'}));
    const rows = div([], 'monomer-rgroups-rows');
    this._form.addElement(rows);
    this._structure.dataset.u2Part = 'error';
    this._structure.hidden = true;
    this.root.append(this._structure, this._form.root);
    const empty = () =>
      span('Draw the monomer with its R groups (R1, R2, ...) to define their caps', 'monomer-rgroups-empty');
    this.runInScope(() => each(rows, this._labels, (label) => `${this._generation}:${label}`,
      (label) => this._row(label), {empty}));
  }

  /** Shows the R groups of a monomer as its library defines them. */
  setDefinitions(rgroups: RGroup[]): void {
    batch(() => {
      this._generation++;
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

    const choice = new ChoiceInput({label, name: label, items: items.peek(), nullable: true,
      tooltipText: `The cap group that closes ${label} when the monomer ends a chain`,
      onChanged: (value) => {
        if (value === OWN && own.peek() !== null)
          def.value = own.peek()!;
        else if (value === null)
          def.value = CapGroups.undefinedRGroup(label);
        else {
          const cap = CapGroups.byName(value);
          if (cap !== null)
            def.value = CapGroups.rgroup(cap, n);
        }
      }});
    // every effect of the row is its controls': they go with the row, not with the editor
    choice.effect(() => {
      const inference = inferred.value;
      if (inference !== null && !inference.known)
        own.value = def.peek();
    });
    choice.effect(() => choice.setItems(items.value));
    choice.follow(() => {
      const inference = inferred.value;
      return inference === null ? null : inference.known ? inference.cap.name : OWN;
    });
    choice.addValidator(() => {
      const problem = this._problemOf(label, def.peek())?.replace(`${label}: `, '') ?? null;
      return problem === null ? null : problem.charAt(0).toUpperCase() + problem.slice(1);
    });
    choice.effect(() => {
      def.value;
      choice.revalidate();
    });

    const info = span('', 'monomer-rgroup-info');
    info.dataset.u2Part = 'description';
    choice.effect(() => {
      const inference = inferred.value;
      const smiles = CapGroups.smilesOf(def.value);
      info.textContent = inference === null ? (smiles || CapGroups.nameOf(def.value) ?
        `not recognized: ${smiles || CapGroups.nameOf(def.value)}` : 'not defined') :
        inference.known ? `${inference.cap.description} · ${CapGroups.numbered(inference.cap.smiles, n)}` +
          (inference.note ? ` (${inference.note})` : '') :
          `${smiles || 'no SMILES'} · ${inference.note ?? 'as defined'}`;
      info.title = info.textContent;
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
        onChanged: (value) => def.value = {...saved.peek(), [FIELD_KEYS[field]]: value.trim()}});
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
    row.dataset.u2Item = label;
    this._form.add(choice, row);
    row.append(fields.root);
    return row;
  }
}
