/* A tumor model the atlas has never seen, planned from its own history (model-history.ts): what a plan keeps of it
   (SavedHistory, inside the plan file so the plan reopens the same way), its registration under its key
   `new:pdx:<id>` so that every stage finds its growth prior by name (priors.ts linePrior), and the notes Plan shows
   about it. The dialog that reads it is history-dialog.ts. */
import {HistoryRow, HistorySource, historyNotes, historyPrior, historyStudy, otherLab, weighHistory} from
  './model-history';
import {HISTORY_PREFIX, LinePrior, enginePriors, pdxHistoryAnchor, setHistoryModel, tumorModels} from './priors';

/** What a plan keeps of a model planned from its history. */
export interface SavedHistory {
  /** the model's readable name, e.g. PDX-LU-117: lists, plans and reports show it (the key stays new:pdx:<id>) */
  name: string;
  source: HistorySource;
  /** the history's control measurements; group means: one row per experiment and day, `animal` the experiment */
  rows: HistoryRow[];
  /** % of the implanted tumors that grew, from the pilot; null: PMX Studio's PDX reference */
  take: number | null;
}

/** What the weighing gave: the growth prior, the weight on the history, its studies and the animals assumed per
 *  group for means. */
export interface HistoryFit {prior: LinePrior, weight: number, studies: number, nAssumed: number}

/** A short stable hash of a text (cyrb53), for the history's key. */
function hashOf(s: string): string {
  let h1 = 0xdeadbeef;
  let h2 = 0x41c6ce57;
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i);
    h1 = Math.imul(h1 ^ c, 2654435761);
    h2 = Math.imul(h2 ^ c, 1597334677);
  }
  h1 = Math.imul(h1 ^ (h1 >>> 16), 2246822507) ^ Math.imul(h2 ^ (h2 >>> 13), 3266489909);
  h2 = Math.imul(h2 ^ (h2 >>> 16), 2246822507) ^ Math.imul(h1 ^ (h1 >>> 13), 3266489909);
  return (4294967296 * (2097151 & h2) + (h1 >>> 0)).toString(36);
}

/** The key the stages know the model by: the same history gives the same key. */
export const historyKey = (h: SavedHistory): string =>
  HISTORY_PREFIX + hashOf(JSON.stringify([h.name, h.source, h.rows, h.take]));

/** The history weighed against the PDX reference, as the engine's growth prior: the PDX lines' mean for what the
 *  history does not measure, their non-progressing share (randomized animals whose tumor stops growing) too. The
 *  pilot's take is not in it: a tumor that never takes is lost before randomization, and the engine reads it from the
 *  saved history (designStudy's preLoss). */
export async function fitHistory(h: SavedHistory): Promise<HistoryFit> {
  const [a, all, models] = await Promise.all([pdxHistoryAnchor(), enginePriors(), tumorModels()]);
  const lines = Object.entries(all.lines).filter(([l]) => models[l]?.class === 'pdx').map(([, p]) => p);
  const nAssumed = a.n_for_group_means;
  const study = historyStudy(h.rows, h.source.meansOnly, h.source.animalsPerGroup, nAssumed);
  if (!study) {
    throw new Error('No growth could be read from the history: each control curve needs at least 4 measurements ' +
      'over which the tumor grows.');
  }
  const w = weighHistory(a, [study], otherLab(h.source));
  return {prior: historyPrior(lines, a, w, study.n), weight: w.weight, studies: 1, nAssumed};
}

const saved = new Map<string, SavedHistory>();
type Fitter = () => Promise<HistoryFit>;
const fits = new Map<string, Fitter>();

/** Registers a model planned from its history under its key (the weighing runs once, when first needed) and
 *  returns the key. */
export function registerHistory(h: SavedHistory): string {
  const key = historyKey(h);
  if (!saved.has(key)) {
    let fit: Promise<HistoryFit> | null = null;
    fits.set(key, () => fit ??= fitHistory(h).catch((e) => {
      fit = null;
      throw e;
    }));
    saved.set(key, h);
    setHistoryModel(key, h.name, () => fits.get(key)!().then((f) => f.prior));
  }
  return key;
}

/** The history of a registered model; null for an atlas model. */
export const savedHistory = (line: string | null | undefined): SavedHistory | null => saved.get(line ?? '') ?? null;

/** The weighing of a registered model's history. */
export const historyFit = (line: string): Promise<HistoryFit> | null => fits.get(line)?.() ?? null;

/** What Plan says about a model planned from its history, where the user decides: where its growth comes from and
 *  how far it was trusted, another lab's wider ranges, the passage and frozen material, assumed animals, a plan
 *  past four weeks, and where the animals to implant come from. None for an atlas model. */
export async function historyNotesOf(line: string, plannedDays: number | null): Promise<string[]> {
  const h = saved.get(line);
  const f = h ? await historyFit(line) : null;
  if (!h || !f)
    return [];
  const take = h.take !== null ? `Animals to implant: from the pilot's take (${h.take}% of the tumors grew): ` +
    `${+(100 - h.take).toFixed(1)}% lost before randomization, or the loss set in Advanced when larger.` :
    'Animals to implant: from the loss before randomization set in Advanced, as the history gives no take.';
  return [...historyNotes(h.source, f.studies, f.weight, f.nAssumed, plannedDays ?? undefined), take];
}

/** Said in Monitor and Analyze for a model planned from its history (D694): their twin takes the PDX lines' twin
 *  fit, which its validation did not cover (it covered Plan's engine). */
export const TWIN_NOTE = 'For this model the animal-by-animal forecast uses how the reference PDX models behave; it ' +
  'has not been checked for this model.';

const fin = (x: unknown): x is number => typeof x === 'number' && Number.isFinite(x);
const answer = (x: unknown) => x === 'yes' || x === 'no' || x === 'unknown';

/** A history read from a plan file someone may have edited: checked before use; null when absent or damaged. */
export function historyOf(x: any): SavedHistory | null {
  const s = x?.source;
  const ok = x && typeof x.name === 'string' && x.name.trim() !== '' && s && answer(s.sameLab) &&
    answer(s.sameHost) && ['early', 'later', 'unknown'].includes(s.passage) && typeof s.frozen === 'boolean' &&
    typeof s.meansOnly === 'boolean' &&
    (s.animalsPerGroup === null || (Number.isInteger(s.animalsPerGroup) && s.animalsPerGroup > 0)) &&
    Array.isArray(x.rows) && x.rows.length > 0 &&
    x.rows.every((r: any) => typeof r?.animal === 'string' && fin(r.day) && fin(r.volume)) &&
    (x.take === null || (fin(x.take) && x.take > 0 && x.take <= 100));
  return !ok ? null : {name: x.name.trim(), source: {sameLab: s.sameLab, sameHost: s.sameHost, passage: s.passage,
    frozen: s.frozen, meansOnly: s.meansOnly, animalsPerGroup: s.animalsPerGroup},
  rows: x.rows.map((r: any) => ({animal: r.animal, day: r.day, volume: r.volume})), take: x.take};
}
