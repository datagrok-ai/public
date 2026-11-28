import {ARCHETYPES, Therapy} from './archetypes';
import {SCHEDULE, scheduleText} from './schedule';

// a text that is a schedule and nothing else
const WHOLE_SCHEDULE = new RegExp(`^\\s*(?:${SCHEDULE.source})\\s*$`, 'i');

export interface Dose {value: number, unit: string}
/** `negative`: the name writes its dose below 0 ("-5 mg/kg"), read as no dose */
export interface ArmPart {drug: string, key: string, dose: Dose | null, schedule: string | null, route: string | null,
  negative?: boolean}

// a dose: a number and a unit (per body weight or area, absolute, radiation, activity, cells)
const DOSE = new RegExp('(\\d+(?:[.,]\\d+)?(?:\\s*[x×]\\s*10\\^?\\d+|e\\d+)?)\\s*(' +
  '(?:mg|g|[uµμ]g|ng)[\\s·]?kg\\^?-1|mg\\/kg|g\\/kg|[uµμ]g\\/kg|ng\\/kg|' +
  'mg\\/m2|mg\\/m²|[nuµμ]mol\\/kg|mmol\\/kg|IU\\/kg|U\\/kg|mg|[uµμ]g|ng|IU|Gy|[kM]Bq|mCi|[uµμ]Ci|PFU|cells|vp)' +
  '(?![a-z])', 'i');
const ROUTE = new RegExp('(?:^|[\\s,;])(i\\.?p\\.?|i\\.?v\\.?|p\\.?o\\.?|s\\.?c\\.?|s\\.?q\\.?|i\\.?t\\.?|' +
  'i\\.?g\\.?|i\\.?m\\.?|oral(?:ly)?|intratumou?ral(?:ly)?|intraperitoneal(?:ly)?|intravenous(?:ly)?|' +
  'subcutaneous(?:ly)?|gavage|peritumou?ral(?:ly)?)(?=$|[\\s,;])', 'i');
const ROUTES: [RegExp, string][] = [[/^i\.?p|intraperit/i, 'i.p.'], [/^i\.?v|intraven/i, 'i.v.'],
  [/^p\.?o|oral|gavage|^i\.?g/i, 'p.o.'], [/^s\.?[cq]|subcut/i, 's.c.'], [/^i\.?t|intratum|peritum/i, 'i.t.'],
  [/^i\.?m/i, 'i.m.']];
// what joins two treatments in a name ("and" / "with" only between spaces; the + of "CD8+ T cells" does not, that
// of "anti-CD8+anti-PD-1" does). In SCHEDULE a repeat count ("x3") is a word of its own, not the X1 of "Cpd-X1".
const JOIN =
  /(?<!\bcd\d{1,3})\s*\+\s*|(?<=\bcd\d{1,3})\+\s*(?=anti|α)|\s+plus\s+|\s+and\s+|\s+with\s+|\s*;\s*|\s+&\s+/i;
// words that are not a drug: the dose's companions
const FILLER = new RegExp('\\b(?:each|both|alone|monotherapy|treatment|combination|combo|at|of|per|doses?|dosed|' +
  'then|supplementation|administration|injection|treated|given)\\b', 'gi');
// what an arm adds as its own control, not a treatment ("Cpd-A + vehicle", "RT + isotype IgG")
const CONTROL =
  /^(?:vehicle|saline|pbs|placebo|control|untreated|dmso|sham)$|isotype|^(?:rat|mouse|hamster|human)?igg\w*$/;
// radiation, however written; also the drug of a part that is only a dose in Gy
const RADIATION = /^(?:rt|radiotherapy|radiation|irradiation|ir|xray|xrt|ionizingradiation)$/;
// where the arm was given, not what: "in Batf3-/- hosts", "on IFNGR1-KO cells"
const CONTEXT = /\s+(?:in|on|into)\s+[^+;]*(?:hosts?|mice|mouse|cells|tumou?rs?|ko|knockout|-\/-|deficient)\b.*$/i;
const SALTS = new RegExp('\\b(?:hydrochloride|hcl|mesylate|mesilate|sodium|potassium|citrate|tartrate|maleate|' +
  'sulfate|sulphate|phosphate|acetate|besylate|tosylate|fumarate|succinate)\\b', 'gi');

/** The controls that turn text right to left or isolate it (U+202E reverses "Cpd-A 5 mg/kg" on screen and in a
 *  document): no name needs them. With `escaped`, also as a JSON text writes them ("‮"). */
const BIDI = /[؜‎‏‪-‮⁦-⁩]/g;
const BIDI_ESCAPED = /[؜‎‏‪-‮⁦-⁩]|\\u(?:061[cC]|200[eEfF]|202[a-eA-E]|206[6-9])/g;
export const noBidi = (s: string, escaped = false) => s.replace(escaped ? BIDI_ESCAPED : BIDI, '');

/** A drug's identity across arms: case, spaces, hyphens, salt forms and the spellings of "anti-" set aside
 *  ("aPD-1", "αPD1", "anti PD-1" are one drug). */
export function drugKey(drug: string): string {
  return drug.toLowerCase().replace(SALTS, '').trim()
    .replace(/^(?:α|a(?=[-\s]?(?:pd|ctla|tim|lag|tigit|cd\d|ox40|4-?1bb|gitr|tgf|vegf|il|csf|nk)))[-\s]?/, 'anti')
    .replace(/^anti[-\s]+/, 'anti').replace(/[\s\-_.,'’]/g, '');
}

/** One treatment's text: the dose, schedule and route found in it, and what is left, the drug. */
function part(text: string): ArmPart | null {
  let rest = ` ${text} `;
  const take = (re: RegExp): RegExpExecArray | null => {
    const m = re.exec(rest);
    if (m)
      rest = rest.slice(0, m.index) + ' ' + rest.slice(m.index + m[0].length);
    return m;
  };
  const d = take(DOSE);
  // a minus written before the dose ("Cpd-A -5 mg/kg", not the hyphen of "Cpd-5"): no dose, and said so (readProblem)
  const negative = !!d && /(?:^|\s)[-−–]\s*$/.test(rest.slice(0, d.index));
  const s = take(SCHEDULE);
  const r = take(ROUTE);
  const route = r ? (ROUTES.find(([re]) => re.test(r[1]))?.[1] ?? r[1]) : null;
  // (a word with no Latin or Greek letter and no digit is a mark, not part of the name: "Trametinib ✓ 中文 ½")
  const drug = rest.replace(FILLER, ' ').replace(/[,;:]+/g, ' ').split(/\s+/)
    .filter((w) => /[A-Za-z0-9Ͱ-Ͽ]/.test(w)).join(' ').trim().replace(/^[-–\s]+|[-–\s]+$/g, '');
  const unitOf = d?.[2] ?? '';
  const named = drug || (/^Gy$/i.test(unitOf) ? 'radiation' : '');
  const key = drugKey(named);
  // (a name with no letter, "1", is no drug)
  if (!/[A-Za-zͰ-Ͽ]/.test(named) || CONTROL.test(key))
    return null;
  const value = d ? Number(d[1].replace(',', '.').replace(/\s*[x×]\s*10\^?(\d+)/, 'e$1')) : NaN;
  return {drug: named, key: RADIATION.test(key) ? 'radiation' : key, schedule: s ? scheduleText(s[0]) : null, route,
    ...negative ? {negative: true} : {},
    dose: d && !negative && Number.isFinite(value) ?
      {value, unit: d[2].replace(/[µμ]/g, 'u').replace(/^(\w+)[\s·]?kg\^?-1$/i, '$1/kg')} : null};
}

/** An arm's treatments, in the order written. A parenthesis is a comment on the treatment before it, not one of its
 *  own ("IL-aC (IL-12 + anti-CTLA-4)" is one treatment); a dose written once for several ("A + B 5 mg/kg each") stays
 *  with the one it is written beside. */
export function parseArmName(name: string): ArmPart[] {
  // a parenthesis that is only a schedule's days ("QD (d0-d4)") is the schedule's, not a comment
  let plain = name.replace(/[([]([^()[\]]*)[)\]]/g,
    (m, inner: string) => WHOLE_SCHEDULE.test(inner) ? ` ${inner} ` : m);
  for (let prev = ''; prev !== plain;)
    [prev, plain] = [plain, plain.replace(/\([^()]*\)|\[[^[\]]*\]/g, ' ')];
  return plain.replace(/\s+/g, ' ').trim().split(JOIN).map((x) => part(x.replace(CONTEXT, '')))
    .filter((x): x is ArmPart => x !== null);
}

/** The control among a study's groups by its name (vehicle, control, PBS, saline, untreated), as Analyze finds it;
 *  undefined: none named so. */
export const controlByName = (groups: string[]): string | undefined =>
  groups.find((g) => /vehicle|control|pbs|saline|untreated/i.test(g));

/** What was understood, in a few words: "Cisplatin · 5 mg/kg · Q7D · i.p. + anti-PD-1 · 10 mg/kg". */
export function describeArm(parts: ArmPart[]): string {
  return parts.map((p) => [p.drug, p.dose && `${p.dose.value} ${p.dose.unit}`, p.schedule, p.route]
    .filter(Boolean).join(' · ')).join(' + ');
}

const ROUTE_WORDS: {[code: string]: string} = {'p.o.': 'by mouth', 'i.p.': 'into the abdomen', 'i.v.': 'into a vein',
  's.c.': 'under the skin', 'i.t.': 'into the tumor', 'i.m.': 'into a muscle'};
const SCHEDULE_WORDS: {[code: string]: string} = {QD: 'once a day', BID: 'twice a day', TID: '3 times a day',
  QOD: 'every other day', QW: 'once a week', BIW: 'twice a week', TIW: '3 times a week', Q7D: 'once a week',
  Q2D: 'every other day'};

/** describeArm's words with its codes in plain words ("Cpd-A · 30 mg/kg · once a day · by mouth"); a code it does
 *  not know stays as written. */
export function plainArm(read: string): string {
  return read.split(' + ').map((p) => p.split(' · ').map((w) => ROUTE_WORDS[w] ?? SCHEDULE_WORDS[w.toUpperCase()] ??
    (/^Q[0-9]+D$/i.test(w) ? `every ${w.slice(1, -1)} days` : w)).join(' · ')).join(' + ');
}

/** What the reading of a treated arm lacks, in a few words ('100 ug: no dose in mg/kg', 'anti-PD-1: no dose' in a
 *  combination); '' when nothing: each drug has a dose per kg, or one that is not a weight (Gy, cells...). */
export function readProblem(parts: ArmPart[]): string {
  if (!parts.length)
    return 'no treatment read in the name';
  return parts.flatMap((p) => {
    if (p.negative)
      return [`${p.drug}: a dose cannot be below 0`];
    if (p.key === 'radiation')
      return [];
    if (!p.dose)
      return [`${parts.length > 1 ? `${p.drug}: ` : ''}dose not given`];
    return /^(?:mg|g|ug|ng)$/i.test(p.dose.unit) ? [`${p.dose.value} ${p.dose.unit}: no dose in mg/kg`] : [];
  }).join('; ');
}

export interface Combination {ab: number, a: number, b: number}

/** Two writings of one drug: the same key, one the start of the other ("antipd1" of "antipd1antibody"), or the same
 *  first word ("siCBX8 knockdown", "caffeine 0.1% drinking water" beside "siCBX8", "caffeine"). */
export function sameDrug(x: ArmPart, y: ArmPart): boolean {
  // the first word, with the next when it is an identifier ("Compound 1" and "Drug B" are not one drug)
  const head = (p: ArmPart) => {
    const w = p.drug.split(' ');
    const id = w.length > 1 && w[1].length <= 3 && /^[A-Za-z]?\d{0,2}[A-Za-z]?$/.test(w[1]);
    return drugKey(id ? `${w[0]} ${w[1]}` : w[0]);
  };
  const [s, l] = x.key.length <= y.key.length ? [x.key, y.key] : [y.key, x.key];
  return x.key === y.key || (s.length >= 3 && l.startsWith(s)) || (head(x).length >= 2 && head(x) === head(y));
}

/** Two doses that do not contradict: one is not written, or both are the same. */
const sameDose = (x: Dose | null, y: Dose | null) => !x || !y ||
  (x.unit.toLowerCase() === y.unit.toLowerCase() && x.value === y.value);

/** Each arm that is the combination of two others: its drugs exactly theirs together, one for one, the two sharing
 *  none (schedules free; a dose that both name must be the same, so with several doses of one drug the arm whose
 *  dose the combination gives is the one found, and none if none does). Arms by index; the first pair found, the arms
 *  read in order. */
export function findCombinations(names: string[]): Combination[] {
  const parts = names.map(parseArmName);
  // the parts of a and b, each paired with a part of ab no other has taken
  const covers = (ab: ArmPart[], ps: ArmPart[]) => {
    const free = [...ab];
    return ps.every((p) => {
      // the combination leaves out a dose of a drug the study gives at another dose too, so which is not said
      const unsaid = (q: ArmPart) => !q.dose && !!p.dose && parts.some((o) => o.length === 1 && o[0] !== p &&
        sameDrug(p, o[0]) && !!o[0].dose && !sameDose(p.dose, o[0].dose));
      const i = free.findIndex((q) => sameDrug(p, q) && sameDose(p.dose, q.dose) && !unsaid(q));
      return i >= 0 && free.splice(i, 1).length > 0;
    });
  };
  const out: Combination[] = [];
  parts.forEach((ab, k) => {
    for (let i = 0; ab.length > 1 && i < parts.length; i++) {
      for (let j = i + 1; j < parts.length; j++) {
        const [a, b] = [parts[i], parts[j]];
        if (i === k || j === k || !a.length || !b.length || a.length + b.length !== ab.length ||
          a.some((p) => b.some((q) => sameDrug(p, q))))
          continue;
        if (covers(ab, [...a, ...b])) {
          out.push({ab: k, a: i, b: j});
          return;
        }
      }
    }
  });
  return out;
}

/** Where a combination `ab` gives a drug otherwise than its arm alone (`monos`' names): its route or schedule, each
 *  written in both ("CPD-B: i.p. q3dx4 here, p.o. QD alone"); [] when none differs. */
export function regimenDifferences(ab: string, monos: string[]): string[] {
  const mine = parseArmName(ab);
  const norm = (s: string | null) => (s ?? '').replace(/[.\s]/g, '').toLowerCase();
  return monos.flatMap(parseArmName).flatMap((p) => {
    const q = mine.find((x) => sameDrug(x, p));
    const keys = (['route', 'schedule'] as const).filter((k) => q && p[k] && q[k] && norm(p[k]) !== norm(q[k]));
    return q && keys.length ?
      [`${q.drug}: ${keys.map((k) => q[k]).join(' ')} here, ${keys.map((k) => p[k]).join(' ')} alone`] : [];
  });
}

/** While an arm's name is typed: a word on combinations when the text points that way (another arm's drug, "+",
 *  "combo", "synergy") and the arm is not yet the combination of two others; else null. */
export function combinationHint(text: string, others: string[], detected: boolean): string | null {
  if (detected || others.length < 2)
    return null;
  const [mine, theirs] = [parseArmName(text), others.flatMap(parseArmName)];
  const pointing = /\+|\bplus\b|\bcombo|\bcombin|\bsynerg/i.test(text) ||
    mine.some((p) => theirs.some((q) => sameDrug(p, q)));
  const two = mine.find((p) => !p.dose && theirs.filter((q) => q.dose && sameDrug(p, q)).length > 1);
  if (pointing && two)
    return `Write the dose of ${two.drug} in the combination's name: the other arms give it at more than one dose.`;
  return pointing ? 'Name it with both drugs, as in "Cpd-A + Cpd-B", and Plan treats it as their combination.' : null;
}

// what kind of treatment a drug's name says, the most specific first. Plan's n depends on it, and reading an
// immunotherapy as a small molecule gives fewer animals than it needs, so an antibody the rules cannot place is read as
// an immunotherapy (the wider spread), the directly acting biologics are the ones named. A name no rule reads is a
// small molecule, and the app says it could not tell; so is a code with only "inhibitor" or "p.o." beside it, which
// says a small molecule but not whether it acts on the tumour or through the immune system. Measured on the atlas's
// curated arms by the preclinical bundle's analysis/therapy_from_name.mjs.
const any = (...rs: RegExp[]) => new RegExp(rs.map((r) => r.source).join('|'));
const IMMUNE_TARGET = any(/pd-?1|pd-?l[12]|ctla-?4|lag-?3|tim-?3|tigit|cd40|ox-?40|4-?1bb|cd137|gitr|cd47|sirp-?a/,
  /csf-?1r|cd73|cd39|cd27|cd28|icos|vista|nkg2a|nkp\d+|cd16|cd3|cd25|cd8|cd4|b7-?h3|il-?\d+[a-z]?r?|tgf-?b\w*|ifn\w*/)
  .source;
const THERAPY_RULES: [Therapy, RegExp][] = [
  ['vaccine_peptide', /vaccin|neoantigen|neo-antigen|peptide|gvax|\bdcs?\b.*pulsed|pulsed.*\bdcs?\b/],
  ['cell_therapy', any(/\bcar[- ]?(?:t|nk|m)\b|\b(?:t|nk|til|cik|cd8\+?|cd4\+?|car)[- ]?cells?\b|adoptive|\btils?\b/,
    /ot-?i+\b|pmel|splenocytes|macrophages?\b|\bcell therapy/)],
  // (a code that starts "si" and a digit reads as an oligo, "STING" only as a word: both lean to more animals)
  ['oligo_rna', /sirna|\bsi-?[a-z]{0,6}\d|mrna|\basos?\b|antisense|oligonucleotide|\bmir-?\d|antagomir|aptamer/],
  ['immune_agonist', any(/\bsting\b|cgamp|\bcdn|adu-?s100|dmxaa|msa-?2|diabzi|\btlr\d*|r848|resiquimod|imiquimod|cpg/,
    /poly ?\(?i:c|polyic|epacadostat|indoximod|\bido\d?\b/)],
  ['biologic', any(/\badcs?\b|vedotin|deruxtecan|\bt-?dxd\b|\bt-?dm1\b|emtansine|govitecan|mmae|ozogamicin/,
    /trastuzumab|cetuximab|bevacizumab|panitumumab|pertuzumab|zanidatamab|ramucirumab|necitumumab|amivantamab|dc101/,
    /(?:anti|α|\ba)-?\s?(?:vegf\w*|her2|erbb\d|egfr|ceacam\d*|trop-?2|cd20|cd19|cd38|cd33|gd2|igf-?1r|met|dll3)\b/,
    /(?:anti|α|\ba)-?\s?(?:nectin-?4|lcn2)\b|fusion protein|\btrap\b|asparaginase|asnase|glucose oxidase/)],
  ['immunotherapy', any(
    new RegExp(String.raw`(?<![a-z0-9])(?:anti-?\s?|alpha-?|α-?|m?pro-?|r[mh]|a-?|m)?(?:${IMMUNE_TARGET})(?![a-z0-9])`),
    /pembrolizumab|nivolumab|ipilimumab|atezolizumab|durvalumab|avelumab|cemiplimab|tremelimumab|tislelizumab/,
    /toripalimab|camrelizumab|sintilimab|relatlimab|magrolimab|abatacept|checkpoint|bispecific|probody/,
    /interleukin|interferon|cytokine|gm-?csf|\btnf/,
    // any other antibody
    /\banti-?\s?[a-z]|\balpha-[a-z]|mab\b|antibod|\big[gea]\d?\b|nanobod/)],
  ['small_molecule', any(/(?:nib|lisib|ciclib|parib|platin|taxel|rubicin|mustine|tecan|citabine|uracil|trexate)\b/,
    /(?:mycin|zomib|sertib|stat|phosphamide)\b|chemo|5-?fu|gemcitabine|fluorouracil|temozolomide|dexamethasone/,
    /metformin|rapamycin|everolimus|tamoxifen/)],
];

/** The more variable of two kinds of treatment (the one that asks for more animals). */
export function moreVariable(a: Therapy, b: Therapy): Therapy {
  return ARCHETYPES[b].kappaSd > ARCHETYPES[a].kappaSd ? b : a;
}

/** The kind of treatment an arm's name says. One drug: its own words, else the whole name (a comment in brackets, a
 *  route). A combination: each drug's own words, and the most variable kind (the cautious n). `placed`: that kind
 *  among the drugs the rules place, null if none; `known`: every drug placed. Not known, `therapy` is `placed` or a
 *  small molecule, and the caller decides the rest (arms-input's therapyFor). Radiation is no drug: beside drugs
 *  their kind is the arm's; alone (`radiation`) it acts on the tumor directly, as a small molecule. */
export function therapyOf(name: string): {therapy: Therapy, known: boolean, placed: Therapy | null,
  radiation: boolean, unplaced: string[]} {
  const all = parseArmName(name);
  const parts = all.filter((p) => p.key !== 'radiation');
  if (all.length && !parts.length)
    return {therapy: 'small_molecule', known: true, placed: 'small_molecule', radiation: true, unplaced: []};
  // a dose in cells is a cell therapy whatever the cells are called
  const own = (p: ArmPart) => kindOf(p.dose?.unit === 'cells' ? `${p.drug} cells` : p.drug);
  const kinds = parts.length > 1 ? parts.map(own) : [(parts[0] && own(parts[0])) ?? kindOf(name)];
  const found = kinds.filter((k): k is Therapy => k !== null);
  const placed = found.length ? found.reduce(moreVariable) : null;
  // the drugs of a combination the name does not place, as written: the card names them
  const unplaced = parts.length > 1 ? parts.filter((_, i) => kinds[i] === null).map((p) => p.drug) : [];
  return {therapy: placed ?? 'small_molecule', known: placed !== null && found.length === kinds.length, placed,
    radiation: false, unplaced};
}

function kindOf(text: string): Therapy | null {
  const t = text.toLowerCase();
  return THERAPY_RULES.find(([, re]) => re.test(t))?.[0] ?? null;
}

/** An arm's name cut after each drug's dose ("Trametinib 1 mg/kg p.o. QD + anti-PD-1 10 mg/kg i.p. BIW" reads
 *  "Trametinib 1 + anti-PD-1 10"); a part with no dose is kept whole. For tight places, the full name on hover. */
export const shortName = (name: string) => name.split(' + ').map((p) => {
  const t = p.trim().split(/\s+/);
  const i = t.findIndex((x, j) => j > 0 && /^\d/.test(x));
  return i < 0 ? p.trim() : t.slice(0, i + 1).join(' ');
}).join(' + ');
