/* The R groups of the monomer in the editor: a chip per R label the structure carries — its cap, red while
   it cannot be saved. A click opens the chip's menu: the common caps, each with what it is and the SMILES
   Save writes, the cap its definition names when that is none of them, and Edit definition, which opens
   the fields a library stores under the chips. */
import {RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';
import {ChoiceChip, Control, Form, TextInput, signal, computed, each, span, div, iconButton, effect, Signal,
  ReadonlySignal, batch, component} from '@datagrok-libraries/u2';
import type {ChoiceChipItem} from '@datagrok-libraries/u2';

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

/** The chip's value for a definition that names a cap outside the common ones. */
const OWN_CAP = '(own)';

export interface RGroupsEditor {
  /** The definitions in label order, as Save writes them: a common cap in its canonical form, any
   * other definition as it stands. */
  readonly rgroups: ReadonlySignal<RGroup[]>;
  /** The first R group that cannot be saved and why, or null. */
  readonly problem: ReadonlySignal<string | null>;
  /** The R group whose definition shows under the chips, or null. */
  readonly open: Signal<string | null>;
  readonly labels: string[];
  /** Shows the R groups of a monomer as its library defines them. */
  setDefinitions(rgroups: RGroup[]): void;
  /** Follows the R labels of the drawn structure: kept labels keep their definitions, a new one
   * starts as the usual cap for its position (R1 H, R2 OH, then H). */
  setLabels(labels: string[]): void;
  /** Says, over the chips, what is wrong with the R groups of the drawn structure; null clears it. */
  setStructureProblem(message: string | null): void;
  /** Sets R`n`'s cap by its name (`OH`) or by its four fields — the way code and copilots edit it. */
  setCap(label: string, cap: string | RGroup): void;
}

export function rgroupsEditor(options: RGroupsEditorOptions = {}): Control & RGroupsEditor {
  const editor = component((self) => {
    const open = signal<string | null>(null);
    const labels = signal<string[]>([]);
    const defs = new Map<string, Signal<RGroup>>();
    /** Bumped by every monomer loaded: its chips are built afresh, while a structure edit keeps them. */
    let generation = 0;
    /** What is wrong with the structure's R groups themselves (an R label twice), over the chips. */
    const structure = span('', 'monomer-rgroups-structure');

    const infer = (def: RGroup): CapInference | null => CapGroups.infer(def, options.canonicalize);

    const saved = (label: string, def: RGroup): RGroup => {
      const inferred = infer(def);
      return inferred?.known ? CapGroups.rgroup(inferred.cap, CapGroups.numberOf(label) ?? 1) :
        {...def, [HELM_RGROUP_FIELDS.LABEL]: label};
    };

    const problemOf = (label: string, def: RGroup): string | null => {
      const inferred = infer(def);
      if (inferred === null)
        return `${label}: pick its cap group, or edit its definition`;
      if (inferred.known)
        return null;
      if (CapGroups.isPlaceholder(inferred.cap))
        return `${label}: X is a placeholder, not a cap group — pick the cap ${label} gets (usually H)`;
      const smiles = CapGroups.smilesOf(def);
      if (smiles === '')
        return `${label}: the cap group ${inferred.cap.name} needs its SMILES`;
      if (options.isValidSmiles && !options.isValidSmiles(smiles))
        return `${label}: the cap SMILES ${smiles} is not valid`;
      return null;
    };

    /** What the chip and its menu say of a definition: the cap, what it is (its friendly name) and the SMILES
     * Save writes. */
    const describe = (label: string, def: RGroup): {cap: string, what: string, smiles: string} => {
      const n = CapGroups.numberOf(label) ?? 1;
      const inference = infer(def);
      const smiles = CapGroups.smilesOf(def);
      if (inference === null) {
        const raw = smiles || CapGroups.nameOf(def);
        return {cap: raw ? '?' : '–', what: raw ? `not recognized: ${raw}` : 'not defined', smiles: ''};
      }
      return inference.known ?
        {cap: inference.cap.name, what: inference.cap.description + (inference.note ? ` (${inference.note})` : ''),
          smiles: CapGroups.numbered(inference.cap.smiles, n)} :
        {cap: inference.cap.name, what: inference.note ?? 'as defined', smiles: smiles || 'no SMILES'};
    };

    /** Why R`n` cannot be saved, without the label the chip already shows. */
    const problemHere = (label: string, def: RGroup): string | null =>
      problemOf(label, def)?.replace(`${label}: `, '') ?? null;

    /** R1 and its cap: a press opens the common caps, each with what it is and the SMILES Save writes, the cap the
     * definition names when it is none of them, and Edit definition. */
    const chipOf = (label: string): HTMLElement => {
      const def = defs.get(label)!;
      const n = CapGroups.numberOf(label) ?? 1;
      // the definition the menu brings back: the latest one that named a cap outside the common ones
      const own = signal<RGroup | null>(null);
      effect(() => {
        const inference = infer(def.value);
        if (inference !== null && !inference.known)
          own.value = def.peek();
      });
      const items = computed((): ChoiceChipItem[] => {
        const caps = COMMON_CAP_GROUPS.map((cap) => ({value: cap.name, label: cap.name,
          description: `${cap.description} · ${CapGroups.numbered(cap.smiles, n)}`}));
        const mine = own.value;
        const ownCap = mine === null ? null : infer(mine)?.cap ?? null;
        if (mine === null || ownCap === null)
          return caps;
        const {what, smiles} = describe(label, mine);
        // a cap of its own that borrows a common name (H over other chemistry) must not read as that one
        return [...caps, {value: OWN_CAP, label: CapGroups.byName(ownCap.name) === null ? ownCap.name :
          `${ownCap.name} (custom)`, description: `${what} · ${smiles}`}];
      });
      const said = computed(() => describe(label, def.value));
      const chip = new ChoiceChip({label, items,
        valueText: computed(() => said.value.cap),
        description: computed(() => said.value.smiles ? `${said.value.what} · ${said.value.smiles}` : said.value.what),
        invalid: computed(() => problemHere(label, def.value)),
        active: computed(() => open.value === label),
        actions: () => {
          const editing = open.peek() === label;
          return [{text: editing ? 'Hide definition' : 'Edit definition…',
            run: () => {
              open.value = editing ? null : label;
            }}];
        },
        onInput: (cap) => {
          def.value = cap === OWN_CAP ? own.peek()! :
            CapGroups.rgroup(COMMON_CAP_GROUPS.find((c) => c.name === cap)!, n);
          options.onInput?.();
        }});
      chip.follow(() => {
        const inference = infer(def.value);
        return inference === null ? null : inference.known ? inference.cap.name : OWN_CAP;
      });
      chip.root.classList.add('monomer-rgroup-chip');
      chip.root.dataset.u2Item = label;
      return chip.root;
    };

    /** The fields a library stores for R`n`, under the chips while its definition is open; they show what
     * Save writes — a cap read from the alternate ID alone shows its SMILES and name too — and an edit
     * starts from that. */
    const definitionOf = (label: string): HTMLElement => {
      const def = defs.get(label)!;
      const written = computed(() => saved(label, def.value));
      const fields = new Form({layout: 'normal'});
      const fieldInput = (caption: string, field: Field, read: (rg: RGroup) => string) => {
        const input = new TextInput({label: caption, name: `${label} ${caption}`, required: true, commitOn: 'change',
          onInput: (value) => {
            def.value = {...written.peek(), [FIELD_KEYS[field]]: value.trim()};
            options.onInput?.();
          }});
        input.follow(() => read(written.value));
        return input;
      };
      fields.addAll([
        fieldInput('Cap SMILES', 'smiles', CapGroups.smilesOf),
        fieldInput('Cap name', 'name', CapGroups.nameOf),
        fieldInput('Alternate ID', 'alternateId', CapGroups.alternateIdOf),
      ]);
      const close = iconButton('times', () => open.value = null, {tooltip: `Close the ${label} definition`});
      const head = div([span(`${label} definition`, 'monomer-rgroup-definition-title'), close],
        'monomer-rgroup-definition-head');
      const box = div([head, fields.root], 'monomer-rgroup-definition');
      box.id = `${editor.rootId}-${label}`;
      box.dataset.u2Item = label;
      fields.effect(() => box.hidden = open.value !== label);
      return box;
    };

    const rgroups = computed(() => labels.value.map((label) => saved(label, defs.get(label)!.value)));
    const problem = computed(() => {
      for (const label of labels.value) {
        const found = problemOf(label, defs.get(label)!.value);
        if (found !== null)
          return found;
      }
      return null;
    });
    structure.dataset.u2Part = 'error';
    structure.hidden = true;
    const chips = div([], 'monomer-rgroups-chips');
    const definitions = div([], 'monomer-rgroups-definitions');
    const root = div([structure, chips, definitions], 'monomer-rgroups-editor');
    const empty = () =>
      span('Draw the monomer with its R groups (R1, R2, ...) to define their caps', 'monomer-rgroups-empty');
    const key = (label: string) => `${generation}:${label}`;
    each(chips, labels, key, (label) => chipOf(label), {empty});
    each(definitions, labels, key, (label) => definitionOf(label));
    // an R group the structure lost takes its definition with it
    self.effect(() => {
      const shown = open.value;
      if (shown !== null && !labels.value.includes(shown))
        open.value = null;
    });

    return {
      root, rgroups, problem, open,
      get labels(): string[] {
        return labels.peek();
      },
      setDefinitions: (monomerRGroups: RGroup[]): void => {
        batch(() => {
          generation++;
          open.value = null;
          defs.clear();
          for (const rg of monomerRGroups) {
            const label = CapGroups.labelOf(rg);
            if (label !== '')
              defs.set(label, signal(rg));
          }
          labels.value = [...defs.keys()].sort(byLabel);
        });
      },
      setLabels: (structureLabels: string[]): void => {
        batch(() => {
          for (const label of structureLabels) {
            if (!defs.has(label))
              defs.set(label, signal(defaultFor(label)));
          }
          labels.value = [...structureLabels].sort(byLabel);
        });
      },
      setStructureProblem: (message: string | null): void => {
        structure.textContent = message ?? '';
        structure.hidden = message === null;
      },
      setCap: (label: string, cap: string | RGroup): void => {
        const def = defs.get(label);
        if (def === undefined)
          throw new Error(`The monomer has no ${label}; its R groups are ${labels.peek().join(', ') || 'none'}`);
        const n = CapGroups.numberOf(label) ?? 1;
        const known = typeof cap === 'string' ? CapGroups.byName(cap) : null;
        if (typeof cap === 'string' && known === null)
          throw new Error(`"${cap}" is not one of the common caps: ${COMMON_CAP_GROUPS.map((c) => c.name).join(', ')}`);
        def.value = known !== null ? CapGroups.rgroup(known, n) :
          typeof cap === 'string' ? CapGroups.undefinedRGroup(label) : {...cap, [HELM_RGROUP_FIELDS.LABEL]: label};
      },
    };
  });
  return editor;
}

function defaultFor(label: string): RGroup {
  const n = CapGroups.numberOf(label) ?? 1;
  const standard = STANDRARD_R_GROUPS.find((rg) => rg.label === label);
  const cap = (standard ? CapGroups.infer(standard)?.cap : null) ?? COMMON_CAP_GROUPS[0];
  return CapGroups.rgroup(cap, n);
}
