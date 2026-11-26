/* Two u2 inputs the Plan form uses in place of the editors generated for two of pmxPlanStudy's parameters:
   - ArmsInput, the `arms` table: one card per treated arm, its name as its title, under it in grey what Plan reads
     in it (the drugs, dose, schedule and route, and the kind of treatment, a click to change it) and the expected
     effect; a remove icon on each and "Add treatment arm" under them. TGI or T/C is the study's (Advanced);
   - CorrectionInput, the `correction` choice: one checkbox (Holm, recommended) with what it means, worded for
     the number of arms declared. Bonferroni stays a value of the function, for scripts and Grokky. */
import * as DG from 'datagrok-api/dg';
import {BoolInput, buttonWithIcon, computed, Control, div, divH, Form, iconButton, Input, InputOptions, Menu,
  NumberInput, ReadonlySignal, Reason, Signal, signal, span, TextArea} from '@datagrok-libraries/u2';
import {THERAPIES, Therapy} from './archetypes';
import {combinationHint, findCombinations, moreVariable, parseArmName, therapyOf} from './arm-name';
import {StudyMeta, detailsName, parseLine} from './study-meta';
import {armColors} from './colors';
import {contextHint, hintLink, noteLine} from './ui';
import {Anchor, FromStudy, SIMPLIFIED, anchorsIn, assumedText, doseText, effectText, fromStudyOf, offerText, rangeAs,
  silence, sizingEffect, undosedText} from './study-effect';

/** `noCombination`: the user said this arm is not the combination the app found in its name.
 *  `therapySource`: "name", Plan reads the kind from the name (and the host); "user", the one chosen (a table
 *  without the column, from a script or an older plan, is the user's). `fromStudy`: the earlier study its effect is
 *  scaled from (`effect` is then the low end of the range, what the arm is sized on) */
export interface Arm {name: string, therapy: Therapy, metric: 'TGI' | 'TC', effect: number,
  therapySource?: 'name' | 'user', fromStudy?: FromStudy, noCombination?: boolean}

/** At most as many arms as the Dunnett tables of crate lsmeans (the sensitivity analysis) and a readable chart hold. */
export const MAX_ARMS = 6;
/** How long typing pauses before the lines under an arm's name read it again (ms). */
const SETTLE_MS = 400;

/** The card of an arm just found to be a combination: its border glows once in `color` and its dot swells, 0.6 s. */
export function combinationPulse(card: HTMLElement, dot: HTMLElement, color: string): void {
  card.dataset.pulsed = 'yes';
  card.animate([{boxShadow: `0 0 0 0 ${color}00`}, {boxShadow: `0 0 0 4px ${color}59`},
    {boxShadow: `0 0 0 0 ${color}00`}], {duration: 600, easing: 'ease-out'});
  dot.animate([{transform: 'scale(1)'}, {transform: 'scale(1.5)'}, {transform: 'scale(1)'}],
    {duration: 600, easing: 'ease-out'});
}

const newArm = (k: number): Arm => ({name: `Treatment ${k}`, therapy: 'small_molecule', metric: 'TGI', effect: 60,
  therapySource: 'name'});

/** A combination's effect when its two drugs act independently (`a`, `b`: their effects, percent, TGI or T/C as
 *  `metric`), rounded: T/C_AB = T/C_A × T/C_B, so TGI_AB = 1 − (1 − TGI_A)(1 − TGI_B) (60 and 60 give 84). The same
 *  line Analyze reads synergy against: Bliss independence on growth rates (synergy.ts, Demidenko & Miller 2019)
 *  adds the two drugs' effects on the rate of log volume, so the logs of their T/C add, and the T/C multiply (Crux
 *  stats.tumor_combination_design's effect_bliss, the kill rates added up, simulated with onset and responders). */
export function independentEffect(a: number, b: number, metric: 'TGI' | 'TC'): number {
  return Math.round(metric === 'TC' ? a * b / 100 : 100 - (100 - a) * (100 - b) / 100);
}

/** What an arm's empty name shows, one after another: what Plan reads in a name (drug, dose, route, schedule) and a
 *  combination. */
export const NAME_EXAMPLES = ['Cpd-A 30 mg/kg QD PO', 'Anti-PD-1 10 mg/kg i.p. Q3D×4',
  'Cpd-A 30 mg/kg QD PO + anti-PD-1 10 mg/kg i.p. Q3D×4'];

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

/** The kinds that act only through the host's own immune system. Not cell therapy: human CAR-T and NK cells are
 *  tested in immunodeficient mice, where nothing rejects them. */
const NEEDS_IMMUNITY: Therapy[] = ['immunotherapy', 'immune_agonist', 'vaccine_peptide'];

/** The drugs of an arm named `name`, of kind `therapy`, that act through an immune system the host `host` lacks
 *  (none on a host that may answer). One drug: by the arm's kind, the user's choice included; a combination: by each
 *  drug's own kind, as its name says it. */
export function lackingImmunity(name: string, therapy: Therapy, host: string | null): string[] {
  if (immuneHost(host))
    return [];
  const parts = parseArmName(name).filter((p) => p.key !== 'radiation');
  if (parts.length > 1) {
    return parts.filter((p) => {
      const k = therapyOf(p.drug).placed;
      return k !== null && NEEDS_IMMUNITY.includes(k);
    }).map((p) => p.drug);
  }
  return NEEDS_IMMUNITY.includes(therapy) ? [parts[0]?.drug ?? (name || 'This arm')] : [];
}

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

/** An arm's name as the study's tables carry it: an arm left unnamed (its card shows the example) is its number. */
export const armName = (a: Arm, i: number) => a.name || `Treatment ${i + 1}`;

/** The arms as the table pmxPlanStudy takes: name, therapy, metric, effect (percent), and where the therapy came
 *  from (therapy_source, read by the form only); with an arm linked to an earlier study, the link as JSON
 *  (from_study, empty for the others); an unticked combination (no_combination). */
export function armsFrame(arms: Arm[]): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('name', arms.map(armName)),
    DG.Column.fromStrings('therapy', arms.map((a) => a.therapy)),
    DG.Column.fromStrings('metric', arms.map((a) => a.metric)),
    DG.Column.fromFloat64Array('effect', Float64Array.from(arms, (a) => a.effect)),
    DG.Column.fromStrings('therapy_source', arms.map((a) => a.therapySource ?? 'user')),
    ...arms.some((a) => a.fromStudy) ?
      [DG.Column.fromStrings('from_study', arms.map((a) => a.fromStudy ? JSON.stringify(a.fromStudy) : ''))] : [],
    // only once someone has unticked a combination: the tables of every other plan stay as they were
    ...(arms.some((a) => a.noCombination) ?
      [DG.Column.fromList('bool', 'no_combination', arms.map((a) => a.noCombination === true))] : []),
  ]);
}

/** A from_study cell's link (none: empty or damaged). */
function linkOf(cell: unknown): FromStudy | undefined {
  try {
    return fromStudyOf(cell ? JSON.parse(String(cell)) : null) ?? undefined;
  } catch {
    return undefined;
  }
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
    ...linkOf(get('from_study', i)) ? {fromStudy: linkOf(get('from_study', i))} : {},
    ...(get('no_combination', i) === true ? {noCombination: true} : {}),
  }));
}

/** Where an arm's effect may come from the user's earlier studies (without it the cards offer none: a script). */
export interface PreviousContext {
  /** the plan's tumor model, and the day its effects are read (null: none yet) */
  line: ReadonlySignal<string | null>;
  evaluationDay: ReadonlySignal<number | null>;
  /** the kind of host mouse a tumor model grows in (its host_class; null: not known) */
  hostOf?: (line: string | null) => string | null;
  /** the earlier study's group that can anchor an arm of this name (previous-studies.ts findAnchor) */
  find: (name: string) => Promise<Anchor | null>;
  /** an earlier study of the arm's drug that says no dose (previous-studies.ts findUndosed) */
  undosed?: (name: string) => Promise<{study: string, group: string} | null>;
  /** the patients' exposure as a mouse dose on the arm's schedule (previous-studies.ts findPatientsDose) */
  patientsDose?: (name: string) => Promise<{text: string, tip: string} | null>;
  /** the study's details (what each group receives), which the name's reading takes first */
  meta?: () => StudyMeta;
  /** an arm's name edited (left, with a new name): its details follow it (study-meta.ts renamedDetails) */
  renamed?: (from: string, to: string) => void;
  /** what Crux makes of it within the drug's PK (study-effect.ts effectFromStudy) */
  estimate: (anchor: Anchor, name: string, day: number) =>
    Promise<Pick<FromStudy, 'point' | 'range' | 'why' | 'anchorTgi' | 'exposure'>>;
  /** changes when tables are opened or closed (a PK study among them): the links are estimated again */
  tables?: ReadonlySignal<number>;
}

/** What the arms take from the study: the tumor model's class, whether effects are TGI or T/C, and the earlier
 *  studies. */
interface StudyContext {host: ReadonlySignal<string | null>, metric: ReadonlySignal<'TGI' | 'TC'>,
  previous?: PreviousContext}

/** What an arm takes from the others: `others`, the names of the arms before it (a new arm that repeats them may be
 *  their combination); `pair`, the two arms it combines, by their drugs, or null; `monos`, those two arms' effects */
interface ArmsContext {others: () => string[], pair: () => string | null, monos: () => [number, number] | null,
  /** the arm's colour in its charts (colors.ts armColors) */
  color: () => string}

/** One arm's card: its fields write through `onChange`; its control is disposed when the arm is removed. */
class ArmCard extends Control {
  /** What holds the arm's fields back, or null; followed by whoever reads it in a signal's scope. */
  problem: () => string | null = () => null;

  constructor(arm: Arm, index: number, study: StudyContext, onChange: (a: Arm) => void, onRemove: () => void,
    ctx: ArmsContext) {
    super();
    this.root.classList.add('pmx-arm');
    Object.assign(this.root.style, {padding: 'var(--dg-space-m)', gap: 'var(--dg-space-s)', display: 'flex',
      flexDirection: 'column'});
    this.build(() => {
      // the name is the card's title, edited in place
      // a box that grows a line at a time, so a long name (a combination with doses and schedules) shows whole;
      // Enter leaves it, as in a one-line field
      const name = new TextArea({name: 'name', value: arm.name, fill: true, autoGrow: true, rows: 1,
        placeholder: `e.g. ${NAME_EXAMPLES[index % NAME_EXAMPLES.length]}`});
      name.root.classList.add('pmx-arm-title');
      const area = name.root.querySelector('textarea')!;
      // what the name is for, without a hover tooltip: one stays up while the name is typed and covers the line under
      // it, the links of what Plan read (the placeholder shows the example, that line what was read)
      area.setAttribute('aria-description', 'The arm\'s name: its drugs, with dose, route and schedule; Plan reads ' +
        `them from it ("${NAME_EXAMPLES[0]}")`);
      // its height is the text's: one line for a short name, no handle to drag
      Object.assign(area.style, {minHeight: '0', resize: 'none', overflow: 'hidden', width: '100%'});
      // the title takes the card's width beside the remove button
      name.root.style.flex = '1 1 auto';
      // the box takes its text's height once laid out and whenever its width changes (a card built before it is
      // shown, a narrower window), not only when typed in
      let width = -1;
      const fit = () => {
        area.style.height = 'auto';
        area.style.height = `${area.scrollHeight}px`;
      };
      // refit on the next frame: resizing inside the observer's own callback is a "ResizeObserver loop" error
      const watch = new ResizeObserver(() => requestAnimationFrame(() => {
        if (area.clientWidth !== width) {
          width = area.clientWidth;
          fit();
        }
      }));
      watch.observe(area);
      this.scope.own(() => watch.disconnect());
      // a name is one line: a pasted line break becomes a space
      area.addEventListener('input', () => {
        if (/[\r\n]/.test(area.value))
          name.value = area.value.replace(/[\r\n]+/g, ' ');
      });
      // the name left changed: this arm's details move from its old name to its new one (the edit says which row)
      let committed = arm.name;
      area.addEventListener('change', () => {
        const now = area.value.replace(/[\r\n]+/g, ' ');
        if (now !== committed && committed)
          study.previous?.renamed?.(committed, now);
        committed = now;
      });
      area.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' && !e.isComposing) {
          e.preventDefault();
          area.blur();
        }
      });
      // what the lines under the name read: the name once typing pauses, so they do not come and go with every key
      // (a half-typed "Trametinib 2 m" reads as another arm and looks for other studies); at once when the name is
      // set without typing (a plan loaded, a test) or the field is left
      const settled = signal(name.peek());
      let pause: ReturnType<typeof setTimeout> | undefined;
      this.effect(() => {
        const n = name.value;
        clearTimeout(pause);
        if (document.activeElement === area)
          pause = setTimeout(() => settled.value = n, SETTLE_MS);
        else
          settled.value = n;
      });
      area.addEventListener('blur', () => {
        clearTimeout(pause);
        settled.value = name.peek();
      });
      this.scope.own(() => clearTimeout(pause));
      name.root.style.fontWeight = '500';
      const remove = iconButton('times', () => onRemove(), {tooltip: 'Remove this arm'});
      const head = divH([name.root, remove]);
      Object.assign(head.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});

      // the kind of treatment: read from the name (and the host) until the user picks one
      const source = signal<'name' | 'user'>(arm.therapySource ?? 'user');
      const picked = signal<Therapy>(arm.therapy);
      const read = computed(() => therapyFor(settled.value, study.host.value));
      const therapy = computed(() => source.value === 'user' ? picked.value : read.value.therapy);
      const pick = (t: Therapy) => {
        picked.value = t;
        source.value = 'user';
      };
      const choose = (e: MouseEvent | KeyboardEvent) => {
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
        // drug · dose · schedule · route, each from the group's details when they say it, else from the name
        const line = parseLine(settled.value, study.previous?.meta?.() ?? null);
        const said: (string | HTMLElement)[] = line ? [line, ' · '] : [];
        if (source.value === 'user')
          return [...said, kindLink(therapy.value, 'Chosen by you: click to change')];
        if (read.value.radiation) {
          return ['Radiation · treated as ', hintLink('acting on the tumor directly', choose,
            'As a small molecule: click to change')];
        }
        if (read.value.known) {
          // a combination is said as one, with its drugs' kinds; it is sized as the more variable of them
          const drugs = parseArmName(settled.value).filter((p) => p.key !== 'radiation');
          if (drugs.length > 1) {
            const kinds = [...new Set(drugs.map((p) => therapyOf(p.drug).placed).filter((k): k is Therapy => !!k))];
            return [...said, hintLink(`combination: ${kinds.map((k) => KIND_WORDS[k]).join(' + ')}`, choose,
              `Read from the name; sized as its more variable drug (${KIND_WORDS[therapy.value]}): click to change`)];
          }
          return [...said, kindLink(therapy.value, 'Read from the name: click to change')];
        }
        if (immuneHost(study.host.value)) {
          // the one question that changes the animals, how a compound the name does not place acts (with the drugs
          // it places, in a combination)
          const way = (label: string, t: Therapy) => {
            const taken = read.value.placed ? moreVariable(read.value.placed, t) : t;
            const l = hintLink(label, () => pick(taken), 'The name does not say. Of the small molecules in the ' +
              'atlas\'s syngeneic studies, about two in three were classed as acting through the immune system, and ' +
              'those need more animals: Plan assumes that until you say');
            // the reading Plan takes now, in bold
            if (therapy.value === taken) {
              const b = document.createElement('b');
              b.append(...l.childNodes);
              l.append(b);
            }
            return l;
          };
          // in a combination, the drug the name does not place named, in the words its own arm's card uses
          const unplaced = read.value.unplaced;
          const what = read.value.placed && unplaced.length ?
            `${unplaced.join(' and ')}: small molecule · works on the tumor ` : 'Small molecule · works on the tumor ';
          return [...said, what, way('directly', 'small_molecule'), ' / ',
            way('through the immune system', 'immune_agonist')];
        }
        return [...said, 'Drug type not recognized: treated as a ',
          kindLink(therapy.value, 'Plan cannot tell the kind of treatment from the name: click to choose it')];
      });
      understood.dataset.pmx = 'understood';
      // a word on combinations while the name is being typed and points there, gone once the field is left; and the
      // combination found, to untick
      const pair = computed(ctx.pair);
      const typing = signal(false);
      name.root.addEventListener('focusin', () => typing.value = true);
      name.root.addEventListener('focusout', () => typing.value = false);
      const hint = contextHint(() => typing.value ?
        combinationHint(settled.value, ctx.others(), pair.value !== null) : null);
      // the tick is the card's own row, a box and its caption in one <label>: the caption in line with the box and a
      // long one wrapping beside it (a u2 BoolInput inside Plan's tall form puts its caption out of line and on one
      // line past the card's edge, and no u2 option changes that)
      const combined = signal(arm.noCombination !== true);
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.name = 'combination';
      this.effect(() => box.checked = combined.value);
      box.addEventListener('change', () => combined.value = box.checked);
      const caption = document.createElement('span');
      this.effect(() => caption.textContent = `Combination of ${pair.value ?? 'two arms'}`);
      const tick = document.createElement('label');
      tick.className = 'pmx-arm-tick';
      // its look here, not in the package CSS (the u2 ratchet caps that): box centred on the caption's first line
      tick.style.cssText = 'display: flex; align-items: flex-start; gap: var(--dg-space-s); line-height: 18px; ' +
        'font-size: var(--dg-font-size); color: var(--dg-text-color); cursor: pointer';
      box.style.cssText = 'flex: none; margin: 2px 0 0; accent-color: var(--dg-accent); cursor: pointer';
      caption.style.cssText = 'min-width: 0; overflow-wrap: anywhere';
      tick.title = 'Untick it if this arm is not given as the combination of those two';
      tick.dataset.pmx = 'combination';
      // a combination found and kept: a thicker border and a dot before the title, in the arm's own colour (its
      // charts'); the tick's words say what it is, the colour draws the eye (inline: the u2 ratchet caps the CSS)
      const dot = document.createElement('span');
      dot.title = 'Detected as a combination of two arms';
      dot.dataset.pmx = 'combination-dot';
      dot.style.cssText = 'flex: none; width: 10px; height: 10px; border-radius: 50%';
      head.prepend(dot);
      // whether it was marked the last time (null: the card's first draw; cards are rebuilt on every arm added,
      // removed or loaded, so only a change on this card, typed or ticked, glows)
      let was: boolean | null = null;
      this.effect(() => {
        const on = pair.value !== null && combined.value;
        const color = ctx.color();
        dot.style.background = color;
        dot.style.display = on ? '' : 'none';
        this.root.dataset.combination = on ? 'yes' : 'no';
        // 1px more border, 1px less padding: the card does not shift
        this.root.style.border = on ? `2px solid ${color}` : '';
        this.root.style.padding = on ? 'calc(var(--dg-space-m) - 1px)' : 'var(--dg-space-m)';
        // found while typing: one brief glow of the border and the dot, once ("this is a combination"); none when
        // the system asks for less motion
        if (was === false && on && !matchMedia('(prefers-reduced-motion: reduce)').matches)
          combinationPulse(this.root, dot, color);
        was = on;
      });
      tick.append(box, caption);
      this.effect(() => tick.style.display = pair.value ? 'flex' : 'none');

      const metricWord = () => study.metric.value === 'TC' ? 'T/C' : 'TGI';
      const effect = new NumberInput({name: 'effect', label: computed(() => `Expected ${metricWord()}`),
        value: arm.effect, min: 5, max: 95, postfix: '%', nullable: false,
        tooltipText: 'The effect expected when the median control reaches the endpoint (TGI or T/C: the study\'s, ' +
          'under Advanced)'});
      // a combination kept: its effect follows its two drugs acting independently until the user types one
      const indep = computed(() => {
        const m = ctx.monos();
        return m && combined.value ? independentEffect(m[0], m[1], study.metric.value) : null;
      });
      const clamped = (v: number) => Math.min(95, Math.max(5, v));
      // ponytail: "typed" is not saved; a loaded combination whose effect is not the independence value counts as
      // typed (a saved plan is never changed behind the user's back), one typed equal to it follows the drugs again
      const typed = signal(indep.peek() !== null && arm.effect !== clamped(indep.peek()!));
      effect.root.addEventListener('input', () => typed.value = true);
      this.effect(() => {
        const v = indep.value;
        if (v !== null && !typed.value && effect.peek() !== clamped(v))
          effect.value = clamped(v);
      });
      const notSynergy = contextHint(() => {
        const [v, x, m] = [indep.value, effect.value, ctx.monos()];
        if (v === null || x === null || !m || !typed.value || (study.metric.value === 'TC' ? x <= v : x >= v))
          return null;
        const line = noteLine(`${study.metric.value === 'TC' ? 'Above' : 'Below'} ${v}%: the two drugs acting ` +
          'independently already give this; not synergy.');
        line.title = `Two drugs acting independently: the combination's T/C is the product of theirs, so ` +
          `${metricWord()} ${m[0]}% and ${m[1]}% give ${v}%. Analyze reads synergy against the same line.`;
        return [line];
      });
      notSynergy.dataset.pmx = 'not-synergy';
      // a drug that works through the immune system, planned in mice that have none: its effect is not to be expected
      const noImmunity = contextHint(() => {
        const drugs = lackingImmunity(settled.value, therapy.value, study.host.value);
        if (!drugs.length)
          return null;
        const [which, it] = [drugs.join(' and '), drugs.length > 1 ? 'them' : 'it'];
        // (a combination's effect from its drugs counts such a drug as inactive: said while it follows them; a group
        // of its own is still sized at the effect typed, so the line says so rather than re-sizing it)
        const inCombination = pair.value && combined.value && !typed.value;
        return [noteLine(`${which} ${drugs.length > 1 ? 'act' : 'acts'} through the immune system, which this ` +
          `model's immunodeficient mice lack: expect little or no effect from ${it} here` +
          (inCombination ? `. The combination's expected effect counts ${it} as inactive.` :
            '; the group is sized as if its expected effect were real.'))];
      });
      noImmunity.dataset.pmx = 'no-immunity';
      // the effect from an earlier study of the same drug: offered quietly under the card, taken with a click
      const prev = study.previous;
      // a name as the arm's line reads it: the plan's details first, then the name
      const asRead = (n: string) => detailsName(n, prev?.meta?.() ?? null);
      // the study's group that can anchor the name, looked up as it changes (the last lookup wins)
      const offer = signal<Anchor | null>(null);
      // or, with none, a study of the drug that says no dose
      const undosed = signal<{study: string, group: string} | null>(null);
      // and where the patients' exposure falls on the mouse's dose scale, from a study that kept both
      const patients = signal<{text: string, tip: string} | null>(null);
      let looked = 0;
      this.effect(() => {
        // (again when the details or the name change)
        const n = asRead(settled.value);
        // (and when the model changes: a study of the plan's own model comes first)
        void prev?.line.value;
        const k = ++looked;
        void prev?.find(n).then((a) => k === looked && (offer.value = a), () => {});
        void prev?.undosed?.(n).then((u) => k === looked && (undosed.value = u), () => {});
        void prev?.patientsDose?.(n).then((d) => k === looked && (patients.value = d), () => {});
      });
      const patientsLine = contextHint(() => {
        const d = patients.value;
        if (!d)
          return null;
        const e = span(d.text);
        e.title = d.tip;
        return [e];
      });
      patientsLine.dataset.pmx = 'patients-dose';
      // the link (as saved, or as taken) and what it gives now: null, the link's own values (a saved plan whose study
      // is not at hand keeps them)
      const link = signal<FromStudy | undefined>(arm.fromStudy);
      const now = signal<Pick<FromStudy, 'point' | 'range' | 'why' | 'doseNew' | 'anchorTgi' | 'exposure'> | null>(
        null);
      const linked = computed(() => link.value && {...link.value, ...now.value ?? {}});
      let asked = 0;
      this.effect(() => {
        const [l, a, n] = [link.value, offer.value, settled.value];
        const [line, day] = [prev?.line.value ?? null, prev?.evaluationDay.value ?? null];
        void prev?.tables?.value;
        const k = ++asked;
        if (!l || !prev)
          return;
        const read = asRead(n);
        const doseNew = parseArmName(read)[0]?.dose?.value ?? l.doseNew;
        // the group it was linked to, even when another dose of the study is now nearer the name's (the offer)
        const mine = a && a.study.name === l.study ?
          a.arm === l.arm ? a : anchorsIn(a.study, read).find((x) => x.arm === l.arm) : undefined;
        if (!mine) {
          now.value = n === arm.name ? null :
            {point: null, range: null, why: 'the name no longer names the drug of your study', doseNew};
          return;
        }
        // (not planned yet: what it gave stands until the plan says its day)
        if (day === null)
          return;
        const why = silence(read, mine, line, day, prev.hostOf);
        if (why) {
          now.value = {point: null, range: null, why, doseNew};
          return;
        }
        // read on the plan's own day, one the study measured (silence)
        void prev.estimate(mine, read, Math.round(day)).then((r) => k === asked && (now.value = {...r, doseNew}));
      });
      const use = () => {
        const a = offer.peek()!;
        now.value = null;
        link.value = {study: a.study.name, arm: a.arm, dosePrev: a.part.dose!.value, unit: a.part.dose!.unit,
          doseNew: parseArmName(asRead(name.peek()))[0]?.dose?.value ?? a.part.dose!.value, point: null, range: null};
      };
      const offered = contextHint(() => {
        const [a, u] = [offer.value, undosed.value];
        // a drug that needs an immune system this model's mice lack: no earlier study speaks for it here (QA: anti-PD-1
        // in a PDX offered its CT26 study, under the warning)
        if (link.value || lackingImmunity(settled.value, therapy.value, study.host.value).length)
          return null;
        if (!a) {
          const drug = parseArmName(asRead(settled.value))[0]?.drug ?? 'this drug';
          return u ? [undosedText(drug, u.study, u.group)] : null;
        }
        return [`${offerText(a, prev?.line.value ?? null)} `, hintLink('Use it to estimate this arm', use,
          `From "${a.study.name}", group ${a.arm}: its effect scaled to this arm's dose`)];
      });
      offered.dataset.pmx = 'offer';
      // the effect a linked arm is sized on, its low end; the field is the user's otherwise
      const sized = computed(() => sizingEffect(linked.value, study.metric.value));
      // how the estimate was made, on hover of its line
      const howTip = (x: FromStudy) => {
        let tip = `Scaled from the effect of ${x.arm}${x.anchorTgi === undefined ? '' :
          ` (TGI ${Math.round(x.anchorTgi)}%)`} in "${x.study}" to ${x.doseNew} ${x.unit}, with the same ` +
          'schedule, route, tumor model and day. The field holds the low end of the range, what the groups are ' +
          'sized on. The rule behind it was checked on 8 dose changes from published studies.';
        // what Plan took from the earlier study because the name does not say it
        const assumed = assumedText(settled.value, x);
        if (assumed)
          tip += ` ${assumed}`;
        // the exposure, as NCA gave it, one line
        if (x.exposure)
          tip += `\n${x.exposure}`;
        // the range's model leaves out the drug's exposure over time (an NCA-derived curve did not pass its check)
        return `${tip}\n${SIMPLIFIED}`;
      };
      // a linked arm: the effect field holds the estimate's low end, what the arm is sized on (typing in it is the
      // user's own value again); one grey line under it says where it comes from, or why there is none
      this.effect(() => {
        const v = sized.value;
        if (v !== null && effect.peek() !== Math.round(v))
          effect.value = Math.round(v);
      });
      effect.root.addEventListener('input', () => {
        if (link.peek())
          link.value = undefined;
      });
      const fromStudy = contextHint(() => {
        const x = linked.value;
        if (!x)
          return null;
        if (!x.range || x.point === null) {
          const said = span(effectText(x, study.metric.value));
          said.title = 'Set the effect by hand.';
          return [said];
        }
        const [lo, hi, point] = rangeAs(x.range, x.point, study.metric.value);
        const said = span(`Scaled from your ${doseText(x)} study and its PK: between ` +
          `${Math.round(lo)} and ${Math.round(hi)}% (≈ ${Math.round(point)}%)`);
        said.title = howTip(x);
        return [said];
      });
      fromStudy.dataset.pmx = 'from-study';

      // a form of its own, captions above: the card keeps its look whatever layout the outer form takes
      const fields = new Form({layout: 'tall'}).addElement(understood).addElement(hint).addElement(tick).add(effect)
        .addElement(fromStudy).addElement(notSynergy).addElement(noImmunity).addElement(offered)
        .addElement(patientsLine);
      // a rename that makes another combination (or none) is not the one that was unticked
      let last = pair.peek();
      this.effect(() => {
        if (pair.value !== last) {
          last = pair.value;
          combined.value = true;
        }
      });
      this.root.append(head, fields.root);
      this.problem = () => effect.validity.value === null ? null : `Arm ${index + 1}: ${effect.validity.value}`;
      // (an empty name stays empty here, so a rebuilt card still shows its example; the tables number it, armName)
      this.effect(() => onChange({name: name.value.trim(), therapy: therapy.value,
        metric: study.metric.value, effect: sized.value ?? effect.value ?? arm.effect, therapySource: source.value,
        ...linked.value ? {fromStudy: linked.value} : {}, ...(combined.value ? {} : {noCombination: true})}));
    });
  }
}

export interface ArmsInputOptions extends InputOptions<DG.DataFrame | null> {
  /** the tumor model's host class (immunocompetent, immunodeficient; null: none chosen or not loaded yet) */
  host?: ReadonlySignal<string | null>;
  /** TGI or T/C, the study's */
  metric?: ReadonlySignal<'TGI' | 'TC'>;
  /** the user's earlier studies, to scale an arm's effect from (none: no offer) */
  previous?: PreviousContext;
}

/** The treated arms, as the `arms` table of pmxPlanStudy. */
export class ArmsInput extends Input<DG.DataFrame | null, ArmsInputOptions> {
  // set in createEditor, which the base constructor calls before this class's field initialisers would run
  private declare _arms: Arm[];
  private declare _cards: Signal<ArmCard[]>;
  private declare _list: HTMLElement;
  private declare _add: HTMLButtonElement;
  // the arms' names as typed, and the arms each name makes the combination of two others
  private declare _names: Signal<string[]>;
  private declare _found: ReadonlySignal<ReturnType<typeof findCombinations>>;
  // the arms' effects, which a combination's follows
  private declare _effects: Signal<number[]>;
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
      const names = this._arms.map(armName);
      const twice = names.find((n, i) => names.indexOf(n) !== i);
      return twice === undefined ? null : `Two arms are named "${twice}": give each arm its own name`;
    });
  }

  protected createEditor(): HTMLElement {
    this._arms = [];
    this._cards = signal<ArmCard[]>([]);
    this._names = signal<string[]>([]);
    this._effects = signal<number[]>([]);
    this._found = computed(() => findCombinations(this._names.value));
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
    // the plan's empty state adds the first arm: this button waits for one
    this.effect(() => this._add.hidden = this.count.value === 0);
    const box = div([this._list, this._add]);
    // its section's title is its caption on screen; a screen reader gets it here
    box.setAttribute('role', 'group');
    box.setAttribute('aria-label', 'Treatment arms');
    // the editor takes the cards' height (u2 gives an editor one control's height)
    Object.assign(box.style, {display: 'flex', flexDirection: 'column', gap: 'var(--dg-space-m)', width: '100%',
      height: 'auto', padding: '0'});
    return box;
  }

  /** Adds a treated arm with an empty name, so its example shows (left empty, it is numbered after the others). */
  add(): void {
    if (this._arms.length < MAX_ARMS)
      this._set([...this._arms, {...newArm(this._arms.length + 1), name: ''}]);
  }

  /** The arms' effects said as `metric` (TGI or T/C): an arm said the other way is converted, T/C = 100 − TGI (the
   *  engine's: T/C_AB = T/C_A × T/C_B against TGI_AB = 1 − (1 − TGI_A)(1 − TGI_B)); one already said so stays. */
  sayAs(metric: 'TGI' | 'TC'): void {
    if (this._arms.some((a) => a.metric !== metric))
      this._set(this._arms.map((a) => a.metric === metric ? a : {...a, metric, effect: 100 - a.effect}));
  }

  /** the arms' names as they stand (each keystroke of a name included) */
  get names(): Signal<string[]> {
    return this._names;
  }

  private _set(arms: Arm[], commit = true): void {
    this._cards.peek().forEach((c) => c.dispose());
    this._arms = arms.slice(0, MAX_ARMS);
    this._names.value = this._arms.map((a) => a.name);
    this._effects.value = this._arms.map((a) => a.effect);
    // the two arms by their drugs ("Cisplatin 5 mg/kg and anti-PD-1"), with the dose when the name has one
    const drugs = (k: number) => parseArmName(this._names.value[k])
      .map((p) => p.dose ? `${p.drug} ${p.dose.value} ${p.dose.unit}` : p.drug).join(' + ');
    const pair = (i: number) => {
      const c = this._found.value.find((x) => x.ab === i);
      return c ? `${drugs(c.a)} and ${drugs(c.b)}` : null;
    };
    // (createEditor runs inside the base constructor: the options are there, the study's signals read from them)
    const study: StudyContext = {host: this.options.host ?? signal(null), metric: this.options.metric ?? signal('TGI'),
      previous: this.options.previous};
    this._cards.value = this._arms.map((a, i) => this.build(() => new ArmCard(a, i, study, (next) => {
      const kind = this._arms[i].therapy;
      this._arms[i] = next;
      this._names.value = this._arms.map((x) => x.name);
      // (a kind changed too: a combination's effect depends on whether its drug needs an immune system)
      if (this._effects.peek()[i] !== next.effect || kind !== next.therapy)
        this._effects.value = this._arms.map((x) => x.effect);
      this._commit();
    },
    () => this._set(this._arms.filter((_, j) => j !== i)),
    {others: () => this._names.value.slice(0, i), pair: () => pair(i), monos: () => {
      const c = this._found.value.find((x) => x.ab === i);
      // a drug that needs an immune system the host lacks counts as inactive in the combination (TGI 0, T/C 100)
      // (the effects read first: a kind changed re-sets them, and this follows)
      const effects = this._effects.value;
      const own = (k: number) => lackingImmunity(this._names.value[k], this._arms[k].therapy, study.host.value)
        .length ? (study.metric.value === 'TC' ? 100 : 0) : effects[k];
      return c ? [own(c.a), own(c.b)] : null;
    }, color: () => armColors(this._names.value, study.previous?.meta?.()?.colors,
      study.previous?.meta?.()?.slots)[this._names.value[i]]})));
    this._list.replaceChildren(...this._cards.peek().map((c) => c.root));
    const full = this._arms.length >= MAX_ARMS ? `At most ${MAX_ARMS} treated arms` : null;
    Reason.apply(this._add, full);
    this._add.title = full ?? 'Add another treated arm';
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
  return a.length === b.length && a.every((x, i) => armName(x, i) === armName(b[i], i) && x.therapy === b[i].therapy &&
    x.metric === b[i].metric && x.effect === b[i].effect &&
    (x.therapySource ?? 'user') === (b[i].therapySource ?? 'user') &&
    JSON.stringify(x.fromStudy ?? null) === JSON.stringify(b[i].fromStudy ?? null) &&
    !x.noCombination === !b[i].noCombination);
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
