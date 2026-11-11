/* Two u2 inputs the Plan form uses in place of the editors generated for two of pmxPlanStudy's parameters:
   - ArmsInput, the `arms` table: one card per treated arm, its name as its title, under it in grey what Plan reads
     in it (the drugs, dose, schedule and route, and the kind of treatment, a click to change it) and the expected
     effect; a remove icon on each and "Add treatment arm" under them. TGI or T/C is the study's (Advanced);
   - CorrectionInput, the `correction` choice: one checkbox (Holm, recommended) with what it means, worded for
     the number of arms declared. Bonferroni stays a value of the function, for scripts and Grokky. */
import * as DG from 'datagrok-api/dg';
import {BoolInput, buttonWithIcon, computed, Control, div, divH, Form, iconButton, Input, InputOptions, Menu,
  NumberInput, ReadonlySignal, Reason, Signal, signal, span, TextInput} from '@datagrok-libraries/u2';
import {THERAPIES, Therapy} from './archetypes';
import {describeArm, moreVariable, parseArmName, therapyOf} from './arm-name';
import {contextHint, hintLink} from './ui';

/** `therapySource`: "name", Plan reads the kind from the name (and the host); "user", the one chosen (a table
 *  without the column, from a script or an older plan, is the user's) */
export interface Arm {name: string, therapy: Therapy, metric: 'TGI' | 'TC', effect: number,
  therapySource?: 'name' | 'user'}

/** At most as many arms as the Dunnett tables of crate lsmeans (the sensitivity analysis) and a readable chart hold. */
export const MAX_ARMS = 6;

const newArm = (k: number): Arm => ({name: `Treatment ${k}`, therapy: 'small_molecule', metric: 'TGI', effect: 60,
  therapySource: 'name'});

/** The kinds of treatment in plain words. */
export const KIND_WORDS: {[t in Therapy]: string} = {
  small_molecule: 'small molecule', biologic: 'antibody or protein acting on the tumor',
  immunotherapy: 'antibody or protein acting through the immune system',
  immune_agonist: 'small molecule acting through the immune system', cell_therapy: 'cell therapy',
  vaccine_peptide: 'vaccine', oligo_rna: 'RNA drug',
};

/** Whether mice of this host class may answer through their immune system: all but immunodeficient ones, and the
 *  unknown host too (no line chosen yet, the atlas still loading), the cautious n. */
export const immuneHost = (host: string | null) => host !== 'immunodeficient';

/** The kind of treatment Plan takes for an arm named `name` on a model whose host class is `host`: the one its name
 *  says. A drug it does not place is, in mice that may answer through the immune system, a small molecule acting
 *  through it (as 211 of the 310 single small molecules curated in the atlas's syngeneic lines, and the cautious n),
 *  else a small molecule; a combination takes the more variable of that and its placed drugs. `known` false: the name
 *  does not place every drug (with an immune host the card then asks how the rest acts). */
export function therapyFor(name: string, host: string | null): ReturnType<typeof therapyOf> {
  const r = therapyOf(name);
  if (r.known || !immuneHost(host))
    return r;
  return {...r, therapy: r.placed ? moreVariable(r.placed, 'immune_agonist') : 'immune_agonist'};
}

/** The arms as the table pmxPlanStudy takes: name, therapy, metric, effect (percent), and where the therapy came
 *  from (therapy_source, read by the form only). */
export function armsFrame(arms: Arm[]): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('name', arms.map((a) => a.name)),
    DG.Column.fromStrings('therapy', arms.map((a) => a.therapy)),
    DG.Column.fromStrings('metric', arms.map((a) => a.metric)),
    DG.Column.fromFloat64Array('effect', Float64Array.from(arms, (a) => a.effect)),
    DG.Column.fromStrings('therapy_source', arms.map((a) => a.therapySource ?? 'user')),
  ]);
}

/** The arms of such a table (none: one default arm, unless `fallback` is off); a missing name is the arm's number.
 *  A value that went through a FuncCall can come back as the platform's own object: DG.toJs gives the DataFrame. */
export function armsOf(value: DG.DataFrame | null | undefined, fallback = true): Arm[] {
  const df: DG.DataFrame | null | undefined = value instanceof DG.DataFrame || value == null ? value : DG.toJs(value);
  if (!df || !(df.rowCount > 0))
    return fallback ? [newArm(1)] : [];
  const get = (c: string, i: number) => df.col(c)?.get(i);
  return Array.from({length: df.rowCount}, (_, i) => ({
    name: String(get('name', i) ?? '').trim() || `Treatment ${i + 1}`,
    therapy: (THERAPIES.includes(get('therapy', i)) ? get('therapy', i) : 'small_molecule') as Therapy,
    metric: get('metric', i) === 'TC' ? 'TC' : 'TGI',
    effect: Number(get('effect', i) ?? 60),
    therapySource: get('therapy_source', i) === 'name' ? 'name' : 'user',
  }));
}

/** What the arms take from the study: the tumor model's class and whether effects are TGI or T/C. */
interface StudyContext {host: ReadonlySignal<string | null>, metric: ReadonlySignal<'TGI' | 'TC'>}

/** One arm's card: its fields write through `onChange`; its control is disposed when the arm is removed. */
class ArmCard extends Control {
  /** What holds the arm's fields back, or null; followed by whoever reads it in a signal's scope. */
  problem: () => string | null = () => null;

  constructor(arm: Arm, index: number, study: StudyContext, onChange: (a: Arm) => void, onRemove: () => void) {
    super();
    this.root.classList.add('pmx-arm');
    Object.assign(this.root.style, {padding: 'var(--dg-space-m)', gap: 'var(--dg-space-s)', display: 'flex',
      flexDirection: 'column'});
    this.build(() => {
      // the name is the card's title, edited in place
      const name = new TextInput({name: 'name', value: arm.name, fill: true, placeholder: `Treatment ${index + 1}`,
        tooltipText: 'The arm\'s name: its drugs, with dose, schedule and route if you like ' +
          '("Cpd-A 10 mg/kg QD p.o.")'});
      name.root.classList.add('pmx-arm-title');
      name.root.style.fontWeight = '500';
      const remove = iconButton('times', () => onRemove(), {tooltip: 'Remove this arm'});
      const head = divH([name.root, remove]);
      Object.assign(head.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});

      // the kind of treatment: read from the name (and the host) until the user picks one
      const source = signal<'name' | 'user'>(arm.therapySource ?? 'user');
      const picked = signal<Therapy>(arm.therapy);
      const read = computed(() => therapyFor(name.value, study.host.value));
      const therapy = computed(() => source.value === 'user' ? picked.value : read.value.therapy);
      const pick = (t: Therapy) => {
        picked.value = t;
        source.value = 'user';
      };
      const choose = (e: MouseEvent) => {
        const menu = Menu.popup();
        for (const t of THERAPIES)
          menu.item(KIND_WORDS[t], () => pick(t), {check: t === therapy.peek()});
        if (source.peek() === 'user')
          menu.separator().item('As read from the name', () => source.value = 'name');
        menu.show({anchor: e.currentTarget as HTMLElement});
      };
      const kindLink = (t: Therapy, tip: string) => hintLink(KIND_WORDS[t], choose, tip);
      // under the name, quietly: what Plan reads in it, and the kind it takes, a click to change it
      const understood = contextHint(() => {
        const parts = parseArmName(name.value);
        const said: (string | HTMLElement)[] = parts.length > 1 || parts.some((p) => p.dose || p.schedule || p.route) ?
          [describeArm(parts), ' · '] : [];
        if (source.value === 'user')
          return [...said, kindLink(therapy.value, 'Chosen by you: click to change')];
        if (read.value.radiation) {
          return ['Radiation · treated as ', hintLink('acting on the tumor directly', choose,
            'As a small molecule: click to change')];
        }
        if (read.value.known)
          return [...said, kindLink(therapy.value, 'Read from the name: click to change')];
        if (immuneHost(study.host.value)) {
          // the one question that changes the animals: how a compound the name does not place acts (with the drugs
          // it places, in a combination)
          const way = (label: string, t: Therapy) => {
            const taken = read.value.placed ? moreVariable(read.value.placed, t) : t;
            const l = hintLink(label, () => pick(taken), 'The name does not say. Most small molecules tested in ' +
              'syngeneic models act through the immune system, and they need more animals: Plan assumes that until ' +
              'you say');
            l.classList.toggle('pmx-hint-current', therapy.value === taken);
            return l;
          };
          const what = read.value.placed ? 'One drug not recognized · it works on the tumor ' :
            'Small molecule · works on the tumor ';
          return [...said, what, way('directly', 'small_molecule'), ' / ',
            way('through the immune system', 'immune_agonist')];
        }
        return [...said, 'Drug type not recognized: treated as a ',
          kindLink(therapy.value, 'Plan cannot tell the kind of treatment from the name: click to choose it')];
      });
      understood.dataset.pmx = 'understood';

      const effect = new NumberInput({name: 'effect',
        label: computed(() => `Expected ${study.metric.value === 'TC' ? 'T/C' : 'TGI'}`),
        value: arm.effect, min: 5, max: 95, postfix: '%', nullable: false,
        tooltipText: 'The effect expected when the median control reaches the endpoint (TGI or T/C: the study\'s, ' +
          'under Advanced)'});
      // a form of its own, captions above: the card keeps its look whatever layout the outer form takes
      const fields = new Form({layout: 'tall'}).addElement(understood).add(effect);
      this.root.append(head, fields.root);
      this.problem = () => effect.validity.value === null ? null : `Arm ${index + 1}: ${effect.validity.value}`;
      this.effect(() => onChange({name: name.value.trim() || `Treatment ${index + 1}`, therapy: therapy.value,
        metric: study.metric.value, effect: effect.value ?? arm.effect, therapySource: source.value}));
    });
  }
}

export interface ArmsInputOptions extends InputOptions<DG.DataFrame | null> {
  /** the tumor model's host class (immunocompetent, immunodeficient; null: none chosen or not loaded yet) */
  host?: ReadonlySignal<string | null>;
  /** TGI or T/C, the study's */
  metric?: ReadonlySignal<'TGI' | 'TC'>;
}

/** The treated arms, as the `arms` table of pmxPlanStudy. */
export class ArmsInput extends Input<DG.DataFrame | null, ArmsInputOptions> {
  // set in createEditor, which the base constructor calls before this class's field initialisers would run
  private declare _arms: Arm[];
  private declare _cards: Signal<ArmCard[]>;
  private declare _list: HTMLElement;
  private declare _add: HTMLButtonElement;
  /** how many arms are declared, for the inputs that word themselves by it */
  declare readonly count: Signal<number>;

  constructor(options: ArmsInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'pmx-arms-input';
    // the study's tables key their rows by arm name: two arms of one name would collapse into one
    this.addValidator(() => {
      const problem = this._cards.value.map((c) => c.problem()).find((p) => p !== null);
      if (problem)
        return problem;
      const names = this._arms.map((a) => a.name);
      const twice = names.find((n, i) => names.indexOf(n) !== i);
      return twice === undefined ? null : `Two arms are named "${twice}"`;
    });
  }

  protected createEditor(): HTMLElement {
    this._arms = [];
    this._cards = signal<ArmCard[]>([]);
    (this as {count: Signal<number>}).count = signal(0);
    this._list = div([], 'pmx-arms');
    Object.assign(this._list.style, {display: 'flex', flexDirection: 'column', gap: 'var(--dg-space-m)'});
    this._add = buttonWithIcon('Add treatment arm', 'plus', () => this.add());
    this._add.dataset.u2 = 'pmx-add-arm';
    // none at first: a new plan starts with no arm (the user adds them)
    this._set(armsOf(this.peek(), false), false);
    // a table written from outside (a run from the history, a script) rebuilds the cards; anything else is
    // u2 funcForm following the call's param back, which reads a table param as its text: ignored
    this.effect(() => {
      const df = this.value as unknown;
      if (df instanceof DG.DataFrame && !sameArms(armsOf(df, false), this._arms))
        this._set(armsOf(df, false), false);
    });
    const box = div([this._list, this._add]);
    // the editor takes the cards' height (u2 gives an editor one control's height)
    Object.assign(box.style, {display: 'flex', flexDirection: 'column', gap: 'var(--dg-space-m)', width: '100%',
      height: 'auto', padding: '0'});
    return box;
  }

  /** Adds a treated arm, numbered after the others. */
  add(): void {
    if (this._arms.length < MAX_ARMS)
      this._set([...this._arms, newArm(this._arms.length + 1)]);
  }

  private _set(arms: Arm[], commit = true): void {
    this._cards.peek().forEach((c) => c.dispose());
    this._arms = arms.slice(0, MAX_ARMS);
    // (createEditor runs inside the base constructor: the options are there, the study's signals read from them)
    const study: StudyContext = {host: this.options.host ?? signal(null), metric: this.options.metric ?? signal('TGI')};
    this._cards.value = this._arms.map((a, i) => this.build(() => new ArmCard(a, i, study, (next) => {
      this._arms[i] = next;
      this._commit();
    },
    () => this._set(this._arms.filter((_, j) => j !== i)))));
    this._list.replaceChildren(...this._cards.peek().map((c) => c.root));
    const full = this._arms.length >= MAX_ARMS;
    Reason.apply(this._add, full ? `At most ${MAX_ARMS} treated arms` : null);
    this._add.title = full ? `At most ${MAX_ARMS} treated arms` : 'Add another treated arm';
    this.count.value = this._arms.length;
    if (commit)
      this._commit();
  }

  private _commit(): void {
    const v = this.peek() as unknown;
    if (!(v instanceof DG.DataFrame) || !sameArms(armsOf(v, false), this._arms))
      this.value = armsFrame(this._arms);
  }
}

function sameArms(a: Arm[], b: Arm[]): boolean {
  return a.length === b.length && a.every((x, i) => x.name === b[i].name && x.therapy === b[i].therapy &&
    x.metric === b[i].metric && x.effect === b[i].effect &&
    (x.therapySource ?? 'user') === (b[i].therapySource ?? 'user'));
}

export interface CorrectionInputOptions extends InputOptions<string | null> {
  /** how many treated arms are declared (ArmsInput.count) */
  arms: ReadonlySignal<number>;
}

/** The multiplicity correction as one checkbox: Holm (recommended) or none, with what it means. */
export class CorrectionInput extends Input<string | null, CorrectionInputOptions> {
  constructor(options: CorrectionInputOptions) {
    super(options, 'holm');
    this.root.dataset.u2 = 'pmx-correction-input';
  }

  protected createEditor(): HTMLElement {
    const box = this.build(() => new BoolInput({label: 'Correct for multiple comparisons', captionAfter: true,
      value: this.peek() !== 'none'}));
    // ponytail: u2's caption-after row tops its caption against a 28px editor, 5px above the box's centre;
    // centred here until u2 centres it
    box.root.style.alignItems = 'center';
    const note = span('', 'pmx-note');
    Object.assign(note.style, {color: 'var(--dg-text-color-muted, #68727f)', fontSize: '11.5px', lineHeight: '1.4'});
    this.build(() => {
      this.effect(() => {
        const v = box.value ? (this.peek() === 'bonferroni' ? 'bonferroni' : 'holm') : 'none';
        if (this.peek() !== v)
          this.value = v;
      });
      this.effect(() => {
        const on = this.value !== 'none';
        if (box.peek() !== on)
          box.value = on;
      });
      // with one treated arm there is nothing to correct: the box is off limits, the reason in its tooltip
      this.effect(() => {
        const k = this.options.arms.value;
        box.disabled = k >= 2 ? null : 'Applies when two or more treated arms are compared with the control';
        note.hidden = k < 2;
        note.textContent = k < 2 ? '' : this.value !== 'none' ?
          `Recommended. ${this.value === 'bonferroni' ? 'Bonferroni' : 'Holm'}'s correction keeps the chance of ` +
          `any false positive over the ${k} comparisons at α; it needs a few more animals.` :
          `Each of the ${k} comparisons is tested at α: fine for an exploratory study, but the chance of at least ` +
          'one false positive grows with the number of arms.';
      });
    });
    const col = div([box.root, note]);
    Object.assign(col.style, {display: 'flex', flexDirection: 'column', gap: '2px', height: 'auto', padding: '0'});
    return col;
  }
}
