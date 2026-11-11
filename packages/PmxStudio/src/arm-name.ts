/* What an arm's free-text name says: its drugs, each with the dose, schedule and route written beside it
   ("Cisplatin 5 mg/kg Q7D i.p. + anti-PD-1 10 mg/kg"), and which arms are the combination of two others.
   Plain functions on strings, no platform: the arms input reads them, and later the PK bridge (dose, schedule). */

import {ARCHETYPES, Therapy} from './archetypes';

export interface Dose {value: number, unit: string}
export interface ArmPart {drug: string, key: string, dose: Dose | null, schedule: string | null, route: string | null}

// a dose: a number and a unit (per body weight or area, absolute, radiation, activity, cells)
const DOSE = new RegExp('(\\d+(?:[.,]\\d+)?(?:\\s*[x×]\\s*10\\^?\\d+|e\\d+)?)\\s*(mg\\/kg|g\\/kg|[uµμ]g\\/kg|ng\\/kg|' +
  'mg\\/m2|mg\\/m²|[nuµμ]mol\\/kg|mmol\\/kg|IU\\/kg|U\\/kg|mg|[uµμ]g|ng|IU|Gy|[kM]Bq|mCi|[uµμ]Ci|PFU|cells|vp)' +
  '(?![a-z])', 'i');
const SCHEDULE = new RegExp('\\b(?:q\\d+(?:d|h|w)(?:\\s*[x×]\\s*\\d+)?|qd|bid|tid|qod|qw|biw|tiw|' +
  '(?:once|twice|three times|\\d+ times) (?:a |per )?(?:day|daily|week|weekly)|daily|weekly|every (?:other|\\d+) ' +
  '(?:day|days|week|weeks)|[x×]\\s*\\d+(?:\\s*doses)?|\\d+\\s*doses)(?![a-z])', 'i');
const ROUTE = new RegExp('(?:^|[\\s,;])(i\\.?p\\.?|i\\.?v\\.?|p\\.?o\\.?|s\\.?c\\.?|s\\.?q\\.?|i\\.?t\\.?|' +
  'i\\.?g\\.?|i\\.?m\\.?|oral(?:ly)?|intratumou?ral(?:ly)?|intraperitoneal(?:ly)?|intravenous(?:ly)?|' +
  'subcutaneous(?:ly)?|gavage|peritumou?ral(?:ly)?)(?=$|[\\s,;])', 'i');
const ROUTES: [RegExp, string][] = [[/^i\.?p|intraperit/i, 'i.p.'], [/^i\.?v|intraven/i, 'i.v.'],
  [/^p\.?o|oral|gavage|^i\.?g/i, 'p.o.'], [/^s\.?[cq]|subcut/i, 's.c.'], [/^i\.?t|intratum|peritum/i, 'i.t.'],
  [/^i\.?m/i, 'i.m.']];
// what joins two treatments in a name ("and" / "with" only between spaces; the + of "CD8+ T cells" does not)
const JOIN = /(?<!\bcd\d{1,3})\s*\+\s*|\s+plus\s+|\s+and\s+|\s+with\s+|\s*;\s*|\s+&\s+/i;
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
  const s = take(SCHEDULE);
  const r = take(ROUTE);
  const route = r ? (ROUTES.find(([re]) => re.test(r[1]))?.[1] ?? r[1]) : null;
  const drug = rest.replace(FILLER, ' ').replace(/[,;:]+/g, ' ').replace(/\s+/g, ' ').trim()
    .replace(/^[-–\s]+|[-–\s]+$/g, '');
  const unitOf = d?.[2] ?? '';
  const named = drug || (/^Gy$/i.test(unitOf) ? 'radiation' : '');
  const key = drugKey(named);
  if (!named || CONTROL.test(key))
    return null;
  const value = d ? Number(d[1].replace(',', '.').replace(/\s*[x×]\s*10\^?(\d+)/, 'e$1')) : NaN;
  return {drug: named, key: RADIATION.test(key) ? 'radiation' : key, schedule: s ? s[0].trim() : null, route,
    dose: d && Number.isFinite(value) ? {value, unit: d[2].replace(/[µμ]/g, 'u')} : null};
}

/** An arm's treatments, in the order written. A parenthesis is a comment on the treatment before it, not one of its
 *  own ("IL-aC (IL-12 + anti-CTLA-4)" is one treatment); a dose written once for several ("A + B 5 mg/kg each") stays
 *  with the one it is written beside. */
export function parseArmName(name: string): ArmPart[] {
  let plain = name;
  for (let prev = ''; prev !== plain;)
    [prev, plain] = [plain, plain.replace(/\([^()]*\)|\[[^[\]]*\]/g, ' ')];
  return plain.replace(/\s+/g, ' ').trim().split(JOIN).map((x) => part(x.replace(CONTEXT, '')))
    .filter((x): x is ArmPart => x !== null);
}

/** What was understood, in a few words: "Cisplatin · 5 mg/kg · Q7D · i.p. | anti-PD-1 · 10 mg/kg". */
export function describeArm(parts: ArmPart[]): string {
  return parts.map((p) => [p.drug, p.dose && `${p.dose.value} ${p.dose.unit}`, p.schedule, p.route]
    .filter(Boolean).join(' · ')).join(' | ');
}

export interface Combination {ab: number, a: number, b: number}

/** Two writings of one drug: the same key, one the start of the other ("antipd1" of "antipd1antibody"), or the same
 *  first word ("siCBX8 knockdown", "caffeine 0.1% drinking water" beside "siCBX8", "caffeine"). */
export function sameDrug(x: ArmPart, y: ArmPart): boolean {
  // the first word, with the next when it is an identifier ("Compound 1", "Drug B": not one drug)
  const head = (p: ArmPart) => {
    const w = p.drug.split(' ');
    const id = w.length > 1 && w[1].length <= 3 && /^[A-Za-z]?\d{0,2}[A-Za-z]?$/.test(w[1]);
    return drugKey(id ? `${w[0]} ${w[1]}` : w[0]);
  };
  const [s, l] = x.key.length <= y.key.length ? [x.key, y.key] : [y.key, x.key];
  return x.key === y.key || (s.length >= 3 && l.startsWith(s)) || (head(x).length >= 2 && head(x) === head(y));
}

/** Each arm that is the combination of two others: its drugs exactly theirs together, one for one, the two sharing
 *  none (doses and schedules free). Arms by index; the first pair found, the arms read in order. */
export function findCombinations(names: string[]): Combination[] {
  const parts = names.map(parseArmName);
  // the parts of a and b, each paired with a part of ab no other has taken
  const covers = (ab: ArmPart[], ps: ArmPart[]) => {
    const free = [...ab];
    return ps.every((p) => {
      const i = free.findIndex((q) => sameDrug(p, q));
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

/** While an arm's name is typed: a word on combinations when the text points that way (another arm's drug, "+",
 *  "combo", "synergy") and the arm is not yet the combination of two others; else null. */
export function combinationHint(text: string, others: string[], detected: boolean): string | null {
  if (detected || others.length < 2)
    return null;
  const theirs = others.flatMap(parseArmName);
  const pointing = /\+|\bplus\b|\bcombo|\bcombin|\bsynerg/i.test(text) ||
    parseArmName(text).some((p) => theirs.some((q) => sameDrug(p, q)));
  return pointing ? 'Name it with both drugs, as in "Cpd-A + Cpd-B", and Plan treats it as their combination.' : null;
}

// what kind of treatment a drug's name says, the most specific first. Plan's n depends on it, and reading an
// immunotherapy as a small molecule gives fewer animals than it needs: an antibody the rules cannot place is read as
// an immunotherapy (the wider spread), the directly acting biologics are the ones named. A name no rule reads is a
// small molecule, and the app says it could not tell. Measured on the atlas's curated arms:
// preclinical bundle analysis/therapy_from_name.mjs.
const IMMUNE_TARGET = /(?:pd-?1|pd-?l[12]|ctla-?4|lag-?3|tim-?3|tigit|cd40|ox-?40|4-?1bb|cd137|gitr|cd47|sirp-?a|csf-?1r|cd73|cd39|cd27|cd28|icos|vista|nkg2a|nkp\d+|cd16|cd3|cd25|cd8|cd4|b7-?h3|il-?\d+[a-z]?r?|tgf-?b\w*|ifn\w*)/.source;
const THERAPY_RULES: [Therapy, RegExp][] = [
  ['vaccine_peptide', /vaccin|neoantigen|neo-antigen|peptide|gvax|\bdcs?\b.*pulsed|pulsed.*\bdcs?\b/],
  ['cell_therapy', /\bcar[- ]?(?:t|nk|m)\b|\b(?:t|nk|til|cik|cd8\+?|cd4\+?|car)[- ]?cells?\b|adoptive|\btils?\b|ot-?i+\b|pmel|splenocytes|macrophages?\b|\bcell therapy/],
  ['oligo_rna', /sirna|\bsi-?[a-z]*\d|mrna|\basos?\b|antisense|oligonucleotide|\bmir-?\d|antagomir|aptamer/],
  ['immune_agonist', /sting|cgamp|\bcdn|adu-?s100|dmxaa|msa-?2|diabzi|\btlr\d*|r848|resiquimod|imiquimod|cpg|poly ?\(?i:c|polyic|epacadostat|indoximod|\bido\d?\b/],
  ['biologic', new RegExp(/\badcs?\b|vedotin|deruxtecan|\bt-?dxd\b|\bt-?dm1\b|emtansine|govitecan|mmae|ozogamicin|/.source +
    /trastuzumab|cetuximab|bevacizumab|panitumumab|pertuzumab|zanidatamab|ramucirumab|necitumumab|amivantamab|dc101|/.source +
    /(?:anti|α|\ba)-?\s?(?:vegf\w*|her2|erbb\d|egfr|ceacam\d*|trop-?2|cd20|cd19|cd38|cd33|gd2|igf-?1r|met|dll3|nectin-?4|lcn2)\b|/.source +
    /fusion protein|\btrap\b|asparaginase|asnase|glucose oxidase/.source)],
  ['immunotherapy', new RegExp(String.raw`(?<![a-z0-9])(?:anti-?\s?|alpha-?|α-?|m?pro-?|r[mh]|a-?|m)?${IMMUNE_TARGET}(?![a-z0-9])|` +
    /pembrolizumab|nivolumab|ipilimumab|atezolizumab|durvalumab|avelumab|cemiplimab|tremelimumab|tislelizumab|/.source +
    /toripalimab|camrelizumab|sintilimab|relatlimab|magrolimab|abatacept|checkpoint|bispecific|probody|/.source +
    /interleukin|interferon|cytokine|gm-?csf|\btnf|/.source +
    // any other antibody
    /\banti-?\s?[a-z]|\balpha-[a-z]|mab\b|antibod|\big[gea]\d?\b|nanobod/.source)],
  ['small_molecule', /(?:nib|lisib|ciclib|parib|platin|taxel|rubicin|mustine|tecan|citabine|uracil|trexate|mycin|zomib|sertib|stat|degrader|protac|phosphamide)\b|inhibitor|agonist|antagonist|chemo|5-?fu|gemcitabine|fluorouracil|temozolomide|dexamethasone|metformin|rapamycin|everolimus|tamoxifen|\bp\.?o\b|oral|gavage|\bi\.?g\b/],
];

/** The kind of treatment an arm's name says, per drug: its own words first, then the whole name (a comment in
 *  brackets, a route); a combination takes its drugs' most variable kind (the cautious n). `known` false: no rule read
 *  it, small molecule by default. */
export function therapyOf(name: string): {therapy: Therapy, known: boolean} {
  const parts = parseArmName(name);
  // a dose in cells is a cell therapy whatever the cells are called
  const own = (p: ArmPart) => kindOf(p.dose?.unit === 'cells' ? `${p.drug} cells` : p.drug);
  const kinds = parts.length > 1 ? parts.map(own) : [(parts[0] && own(parts[0])) ?? kindOf(name)];
  const known = kinds.filter((k): k is Therapy => k !== null);
  if (!known.length)
    return {therapy: 'small_molecule', known: false};
  return {therapy: known.reduce((a, b) => ARCHETYPES[b].kappaSd > ARCHETYPES[a].kappaSd ? b : a), known: true};
}

function kindOf(text: string): Therapy | null {
  const t = text.toLowerCase();
  return THERAPY_RULES.find(([, re]) => re.test(t))?.[0] ?? null;
}
