/* Two u2 inputs the Plan form uses in place of the editors generated for two of pmxPlanStudy's parameters:
   - ArmsInput, the `arms` table: one card per treated arm (name, therapy, metric, expected effect), a remove icon
     on each and "Add treatment arm" under them;
   - CorrectionInput, the `correction` choice: one checkbox (Dunnett, recommended) with what it means, worded for
     the number of arms declared. Bonferroni stays a value of the function, for scripts and Grokky. */
import * as DG from 'datagrok-api/dg';
import {BoolInput, buttonWithIcon, ChoiceInput, Control, div, divH, Form, iconButton, Input, InputOptions, NumberInput,
  RadioInput, ReadonlySignal, Signal, signal, span, TextInput} from '@datagrok-libraries/u2';
import {ARCHETYPES, THERAPIES, Therapy} from './archetypes';

export interface Arm {name: string, therapy: Therapy, metric: 'TGI' | 'TC', effect: number}

/** At most as many arms as the Dunnett tables of crate lsmeans and a readable chart hold. */
export const MAX_ARMS = 6;

const newArm = (k: number): Arm => ({name: `Treatment ${k}`, therapy: 'small_molecule', metric: 'TGI', effect: 60});

/** The arms as the table pmxPlanStudy takes: name, therapy, metric, effect (percent). */
export function armsFrame(arms: Arm[]): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('name', arms.map((a) => a.name)),
    DG.Column.fromStrings('therapy', arms.map((a) => a.therapy)),
    DG.Column.fromStrings('metric', arms.map((a) => a.metric)),
    DG.Column.fromFloat64Array('effect', Float64Array.from(arms, (a) => a.effect)),
  ]);
}

/** The arms of such a table (none: one default arm); a missing name is the arm's number. A value that went
 *  through a FuncCall can come back as the platform's own object: DG.toJs gives the DataFrame. */
export function armsOf(value: DG.DataFrame | null | undefined): Arm[] {
  const df: DG.DataFrame | null | undefined = value instanceof DG.DataFrame || value == null ? value : DG.toJs(value);
  if (!df || !(df.rowCount > 0))
    return [newArm(1)];
  const get = (c: string, i: number) => df.col(c)?.get(i);
  return Array.from({length: df.rowCount}, (_, i) => ({
    name: String(get('name', i) ?? '').trim() || `Treatment ${i + 1}`,
    therapy: (THERAPIES.includes(get('therapy', i)) ? get('therapy', i) : 'small_molecule') as Therapy,
    metric: get('metric', i) === 'TC' ? 'TC' : 'TGI',
    effect: Number(get('effect', i) ?? 60),
  }));
}

/** One arm's card: its fields write through `onChange`; its control is disposed when the arm is removed. */
class ArmCard extends Control {
  constructor(arm: Arm, index: number, removable: boolean, onChange: (a: Arm) => void, onRemove: () => void) {
    super();
    this.root.classList.add('u2-card', 'pmx-arm');
    Object.assign(this.root.style, {padding: 'var(--dg-space-m)', gap: 'var(--dg-space-s)', display: 'flex',
      flexDirection: 'column'});
    this.runInScope(() => {
      const name = new TextInput({label: 'Name', value: arm.name, fill: true});
      const therapy = new ChoiceInput({label: 'Therapy', value: arm.therapy, fill: true, nullable: false,
        items: THERAPIES.map((t) => ({value: t, label: ARCHETYPES[t].label})),
        tooltipText: 'The kind of treatment this arm gives'});
      // T/C on screen, TC in the table the function reads
      const metric = new RadioInput({label: 'Effect as', items: ['TGI', 'T/C'], buttons: true,
        value: arm.metric === 'TC' ? 'T/C' : 'TGI',
        itemTooltips: {'TGI': 'Tumor growth inhibition', 'T/C': 'Treated / control volume ratio'}});
      const effect = new NumberInput({label: 'Expected', value: arm.effect, min: 5, max: 95, postfix: '%',
        tooltipText: 'The effect expected when the median control reaches the endpoint'});
      const remove = iconButton('times', () => onRemove(), {tooltip: removable ? 'Remove this arm' :
        'A study needs at least one treated arm'});
      remove.disabled = !removable;
      const head = divH([span(`Arm ${index + 1}`, 'pmx-arm-title'), remove]);
      Object.assign(head.style, {justifyContent: 'space-between', alignItems: 'center', fontWeight: '500'});
      // a form of its own, captions above: the card keeps its look whatever layout the outer form takes
      const fields = new Form({layout: 'tall'}).add(name).add(therapy).add(metric).add(effect);
      this.root.append(head, fields.root);
      this.effect(() => onChange({name: name.value.value.trim() || `Treatment ${index + 1}`,
        therapy: (therapy.value.value ?? 'small_molecule') as Therapy,
        metric: metric.value.value === 'T/C' ? 'TC' : 'TGI', effect: effect.value.value ?? arm.effect}));
    });
  }
}

/** The treated arms, as the `arms` table of pmxPlanStudy. */
export class ArmsInput extends Input<DG.DataFrame | null> {
  // set in createEditor, which the base constructor calls before this class's field initialisers would run
  private declare _arms: Arm[];
  private declare _cards: ArmCard[];
  private declare _list: HTMLElement;
  private declare _add: HTMLButtonElement;
  /** how many arms are declared, for the inputs that word themselves by it */
  declare readonly count: Signal<number>;

  constructor(options: InputOptions<DG.DataFrame | null>) {
    super(options, null);
    this.root.dataset.u2 = 'pmx-arms-input';
  }

  protected createEditor(): HTMLElement {
    this._arms = [];
    this._cards = [];
    (this as {count: Signal<number>}).count = signal(1);
    this._list = div([], 'pmx-arms');
    Object.assign(this._list.style, {display: 'flex', flexDirection: 'column', gap: 'var(--dg-space-m)'});
    this._add = buttonWithIcon('Add treatment arm', 'plus',
      () => this._set([...this._arms, newArm(this._arms.length + 1)]));
    this._add.dataset.u2 = 'pmx-add-arm';
    this._set(armsOf(this.value.peek()), false);
    // a table written from outside (a run from the history, a script) rebuilds the cards; anything else is
    // u2 funcForm following the call's param back, which reads a table param as its text: ignored
    this.effect(() => {
      const df = this.value.value as unknown;
      if (df instanceof DG.DataFrame && !sameArms(armsOf(df), this._arms))
        this._set(armsOf(df), false);
    });
    const box = div([this._list, this._add]);
    // the editor takes the cards' height (u2 gives an editor one control's height)
    Object.assign(box.style, {display: 'flex', flexDirection: 'column', gap: 'var(--dg-space-m)', width: '100%',
      height: 'auto', padding: '0'});
    return box;
  }

  private _set(arms: Arm[], commit = true): void {
    this._cards.forEach((c) => c.dispose());
    this._arms = arms.slice(0, MAX_ARMS);
    this._cards = this._arms.map((a, i) => this.runInScope(() => new ArmCard(a, i, this._arms.length > 1,
      (next) => {
        this._arms[i] = next;
        this._commit();
      },
      () => this._set(this._arms.filter((_, j) => j !== i)))));
    this._list.replaceChildren(...this._cards.map((c) => c.root));
    this._add.disabled = this._arms.length >= MAX_ARMS;
    this._add.title = this._add.disabled ? `At most ${MAX_ARMS} treated arms` : 'Add another treated arm';
    this.count.value = this._arms.length;
    if (commit)
      this._commit();
  }

  private _commit(): void {
    const v = this.value.peek() as unknown;
    if (!(v instanceof DG.DataFrame) || !sameArms(armsOf(v), this._arms))
      this.value.value = armsFrame(this._arms);
  }
}

function sameArms(a: Arm[], b: Arm[]): boolean {
  return a.length === b.length && a.every((x, i) => x.name === b[i].name && x.therapy === b[i].therapy &&
    x.metric === b[i].metric && x.effect === b[i].effect);
}

export interface CorrectionInputOptions extends InputOptions<string | null> {
  /** how many treated arms are declared (ArmsInput.count) */
  arms: ReadonlySignal<number>;
}

/** The multiplicity correction as one checkbox: Dunnett (recommended) or none, with what it means. */
export class CorrectionInput extends Input<string | null, CorrectionInputOptions> {
  constructor(options: CorrectionInputOptions) {
    super(options, 'dunnett');
    this.root.dataset.u2 = 'pmx-correction-input';
  }

  protected createEditor(): HTMLElement {
    const box = new BoolInput({label: 'Correct for multiple comparisons', captionAfter: true,
      value: this.value.peek() !== 'none'});
    // ponytail: u2's caption-after row tops its caption against a 28px editor, 5px above the box's centre;
    // centred here until u2 centres it
    box.root.style.alignItems = 'center';
    const note = span('', 'pmx-note');
    Object.assign(note.style, {color: 'var(--dg-text-color-muted, #68727f)', fontSize: '11.5px', lineHeight: '1.4'});
    this.runInScope(() => {
      this.effect(() => {
        const v = box.value.value ? (this.value.peek() === 'bonferroni' ? 'bonferroni' : 'dunnett') : 'none';
        if (this.value.peek() !== v)
          this.value.value = v;
      });
      this.effect(() => {
        const on = this.value.value !== 'none';
        if (box.value.peek() !== on)
          box.value.value = on;
      });
      // with one treated arm there is nothing to correct: the box is off limits, the reason in its tooltip
      this.effect(() => {
        const k = this.options.arms.value;
        box.enabled = k >= 2;
        box.tooltipText = k < 2 ? 'Applies when two or more treated arms are compared with the control' : '';
        note.hidden = k < 2;
        note.textContent = k < 2 ? '' : this.value.value !== 'none' ?
          `Recommended. Dunnett's test keeps the chance of any false positive over the ${k} comparisons at α; it ` +
          'needs a few more animals.' :
          `Each of the ${k} comparisons is tested at α: fine for an exploratory study, but the chance of at least ` +
          'one false positive grows with the number of arms.';
      });
    });
    const col = div([box.root, note]);
    Object.assign(col.style, {display: 'flex', flexDirection: 'column', gap: '2px', height: 'auto', padding: '0'});
    return col;
  }
}
