import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

import {DesirabilityProfile, MpoCalculator, desirabilityKey} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {BaseWeights, campoWeights, explainedShare, finalWeights, pearsonMatrix, prepareDesirabilities,
  robustInverse} from './campo-math';
import {CampoViewers, WeightStatus, weightsFrame} from './campo-viewers';
import {applyDesirabilityTags} from './utils';

const DEBOUNCE_MS = 150;

function listOf(names: string[], max = 4): string {
  return names.length <= max ? names.join(', ') : `${names.slice(0, max).join(', ')} and ${names.length - max} more`;
}

export type CampoMessageKind = 'error' | 'warning';

export interface CampoHost {
  readonly df: DG.DataFrame | null;
  readonly profile: DesirabilityProfile;
  readonly editor: MpoProfileEditor;
  readonly calculator: MpoCalculator | null;
  showCampoMessage(message: string | null, kind: CampoMessageKind): void;
}

/** Everything stage A depends on; weights, priorities and pins are not part of it. */
export function stageAKey(df: DG.DataFrame, profile: DesirabilityProfile,
  mapping: Record<string, string | null>): string {
  return df.id + '|' + Object.entries(profile.properties).map(([name, p]) => {
    const col = mapping[name] ?? name;
    return `${name}>${col}>${df.col(col) ? 1 : 0}>${desirabilityKey(p)}`;
  }).join('|');
}

/** The last stage A result, kept for the "why this weight" tooltip and the AI state. */
interface StageA {
  names: string[];
  R: Float64Array;
  Rinv: Float64Array;
  active: Int32Array;
  weights: BaseWeights;
  /** Properties with no column in the dataset. */
  unmapped: Set<number>;
  /** Property index → index of the earlier property that already uses its column. */
  sharedWith: Map<number, number>;
  columnOf: (string | null)[];
}

export interface CampoState {
  raw: Record<string, number | null>;
  final: Record<string, number>;
  status: Record<string, WeightStatus>;
  /** Why the latest priorities and pins were not applied; the weights above are the last ones that were. */
  error: string | null;
}

/** Correlation-aware weights for the profile view: stage A (desirabilities → R → R⁻¹ → base weights) reruns when
 *  the data or the curves change, stage B (priorities, pins) when the user edits those. */
export class CampoController {
  private stageA: StageA | null = null;
  private final: Float64Array | null = null;
  private status: WeightStatus[] = [];
  private error: string | null = null;
  private keyA = '';
  private warningsA: string[] = [];
  private subs: Subscription[] = [];

  /** `viewers` is null without a dataset. */
  constructor(private host: CampoHost, private viewers: CampoViewers | null) {
    const editor = host.editor;
    // curve dragging fires many events
    this.subs.push(editor.onChanged.pipe(debounceTime(DEBOUNCE_MS)).subscribe(() => this.runStageA()));
    this.subs.push(editor.onWeightingInputChanged.subscribe(() => this.update()));
    if (host.df)
      this.subs.push(host.df.onColumnsAdded.pipe(debounceTime(DEBOUNCE_MS)).subscribe(() => this.runStageA()));
    editor.explainWeight = (name) => this.why(name);
  }

  private get names(): string[] {
    return Object.keys(this.host.profile.properties);
  }

  /** Brings the weights up to date now: stage A if what it depends on changed since its last run (a rename or a
   *  deletion still inside the debounce), stage B otherwise. */
  update(): void {
    const key = this.keyA;
    this.runStageA();
    if (this.keyA === key)
      this.runStageB();
  }

  private fail(message: string): void {
    this.stageA = null;
    this.final = null;
    this.warningsA = [];
    this.viewers?.hide();
    this.host.showCampoMessage(message, 'error');
  }

  runStageA(force = false): void {
    const {df, profile, editor} = this.host;
    const names = this.names;
    if (names.length < 2) {
      this.keyA = '';
      editor.setWeights(Object.fromEntries(names.map((nm) => [nm, 1])), null);
      return this.fail('Needs at least two properties');
    }
    if (!df) {
      this.keyA = '';
      return this.fail('Select a dataset to compute weights');
    }
    const key = stageAKey(df, profile, editor.columnMapping);
    if (!force && key === this.keyA)
      return;
    this.keyA = key;

    // a column already used by an earlier property is not used again: one desirability per column
    const seen = new Map<DG.Column, number>();
    const unmapped = new Set<number>();
    const sharedWith = new Map<number, number>();
    const columns = names.map((nm, k) => {
      const c = df.col(editor.columnMapping[nm] ?? nm);
      if (!c)
        unmapped.add(k);
      else if (seen.has(c))
        sharedWith.set(k, seen.get(c)!);
      else {
        seen.set(c, k);
        return c;
      }
      return null;
    });
    const tagged = applyDesirabilityTags(df, Object.fromEntries(names
      .map((nm, k) => [columns[k]?.name, profile.properties[nm]]).filter(([c]) => c)), true);
    const maps = (this.host.calculator ?? new MpoCalculator()).desirabilities(df, tagged);
    const D = columns.map((c) => c ? maps[tagged.indexOf(c)].D : null);

    const prep = prepareDesirabilities(D, df.rowCount);
    const m = prep.active.length;
    if (prep.rows === 0)
      return this.fail('No rows with all properties present');
    if (m < 2 && prep.rows > names.length) {
      const flat = names.filter((_, k) => columns[k] && !prep.active.includes(k));
      if (flat.length)
        return this.fail(`Desirability does not vary for ${listOf(flat)}: shape their curves to compute weights`);
    }
    if (m < 2 || prep.rows <= m)
      return this.fail('Not enough data to compute correlation-aware weights');

    const activeNames = Array.from(prep.active, (k) => names[k]);
    this.viewers?.showCorr(DG.DataFrame.fromColumns(prep.cols.map((c, j) =>
      DG.Column.fromFloat32Array(activeNames[j], c))), activeNames);
    const R = pearsonMatrix(prep.cols, prep.rows);
    const inv = robustInverse(R, m);
    if (!inv)
      return this.fail('Correlation matrix is singular');

    this.warningsA = [];
    if (prep.rows < 10 * m)
      this.warningsA.push(`Only ${prep.rows} complete rows for ${m} properties: weights may be unreliable.`);
    if (prep.rows < df.rowCount / 2)
      this.warningsA.push(`Weights are derived from ${prep.rows} of ${df.rowCount} rows that have every property.`);
    if (inv.lambda > 0)
      this.warningsA.push(`Correlation matrix is nearly singular: regularised with λ = ${inv.lambda}.`);
    this.stageA = {names, R, Rinv: inv.Rinv, active: prep.active, unmapped, sharedWith,
      weights: campoWeights(inv.Rinv, m, prep.active, names.length), columnOf: columns.map((c) => c?.name ?? null)};
    this.runStageB();
  }

  runStageB(): void {
    if (!this.stageA)
      return;
    const {names, weights} = this.stageA;
    const props = names.map((nm) => this.host.profile.properties[nm]);
    const prio = Float64Array.from(props, (p) => p.priority ?? 1);
    const pinned = Float64Array.from(props, (p) => p.pinnedWeight ?? NaN);
    const {final, error} = finalWeights(weights.base, prio, pinned);
    this.error = error;
    if (!final)
      return this.host.showCampoMessage(error, 'error');
    this.host.showCampoMessage(this.warningsA.length ? this.warningsA.join(' ') : null, 'warning');

    this.final = final;
    this.status = names.map((_, i) => !Number.isNaN(pinned[i]) ? 'pinned' : final[i] === 0 ? 'dropped' : 'kept');
    this.viewers?.showWeights(weightsFrame(names, weights.raw, final, this.status), weights.raw, final,
      weights.clipped);
    this.host.editor.setWeights(Object.fromEntries(names.map((nm, i) => [nm, final[i]])),
      Object.fromEntries(names.map((nm, i) => [nm, weights.base[i]])));
  }

  /** Raw and final weights by property, for the AI state; null until stage B has run. */
  get state(): CampoState | null {
    if (!this.stageA || !this.final)
      return null;
    const {names, weights} = this.stageA;
    const byName = <T>(f: (i: number) => T) => Object.fromEntries(names.map((nm, i) => [nm, f(i)]));
    return {
      raw: byName((i) => Number.isNaN(weights.raw[i]) ? null : weights.raw[i]),
      final: byName((i) => this.final![i]),
      status: byName((i) => this.status[i]),
      error: this.error,
    };
  }

  /** The "why this weight" tooltip of a property name. */
  private why(name: string): HTMLElement | null {
    const a = this.stageA;
    const i = a?.names.indexOf(name) ?? -1;
    if (!a || i < 0 || !this.final)
      return null;
    const lines = [ui.divText(`${name} · final weight ${this.final[i].toFixed(3)}`, 'chem-mpo-why-title')];
    const pos = a.active.indexOf(i);
    if (pos < 0) {
      const shared = a.sharedWith.get(i);
      lines.push(ui.divText(a.unmapped.has(i) ? 'No column in this dataset — excluded.' :
        shared != null ? `Column ${a.columnOf[shared]} is already used by ${a.names[shared]} — excluded.` :
          'Constant in this dataset — excluded.'));
      return ui.divV(lines);
    }

    const m = a.active.length;
    const r2 = explainedShare(a.Rinv, m, pos);
    lines.push(ui.divText(`${Math.round(r2 * 100)}% of its desirability is explained by the other properties.`));
    const pinned = this.host.profile.properties[name]?.pinnedWeight;
    if (pinned != null)
      lines.push(ui.divText(`Pinned at ${pinned.toFixed(3)}.`));
    else if (a.weights.raw[i] < 0)
      lines.push(ui.divText(`Raw weight ${a.weights.raw[i].toFixed(2)} < 0, so it was set to 0.`));

    const partners = Array.from(a.active, (k, q) => ({k, r: a.R[pos * m + q]}))
      .filter((o) => o.k !== i)
      .sort((x, y) => Math.abs(y.r) - Math.abs(x.r))
      .slice(0, 3);
    lines.push(ui.divText('Most correlated', 'chem-mpo-why-caption'));
    for (const o of partners) {
      const bar = ui.div([], `chem-mpo-why-bar ${o.r >= 0 ? 'chem-mpo-why-bar-pos' : 'chem-mpo-why-bar-neg'}`);
      bar.style.width = `${Math.round(Math.abs(o.r) * 40)}px`;
      lines.push(ui.divH([ui.divText(a.names[o.k], 'chem-mpo-why-name'), bar, ui.divText(o.r.toFixed(2))],
        'chem-mpo-why-row'));
    }
    return ui.divV(lines);
  }

  dispose(): void {
    for (const s of this.subs)
      s.unsubscribe();
    this.subs = [];
    this.host.editor.explainWeight = null;
    this.stageA = null;
    this.final = null;
    this.keyA = '';
  }
}
