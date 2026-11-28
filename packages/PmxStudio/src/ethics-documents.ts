/* Plan's documents for an ethics application, as Word files the user completes: the study plan (laid out after the
   ARRIVE study plan, NC3Rs) and the non-technical summary (EU 2020/569, Annex I Part A, labels in English or
   Spanish). Each is a form of fields: what the plan says, or a grey "[…]" for what only the user knows. The same
   form gives the checklist dialog's rows (one per block) and the document. They take the plan's facts, not where they
   are stored. */
import {ArriveItem, receivesLine} from './arrive';
import {PLAN_WEIGHT_LIMIT} from './ethics-text';
import {Block, Cell, Run, fill} from './docx';
import {StudyMeta, animalsLine, dateOf} from './study-meta';

/** A field: its label and what the plan says (null: the user fills it). */
export type Field = [label: string, value: Cell | null];
/** A block of fields: one row of the dialog; `fix` where the app fills the gap (the study's details). */
/** `short`: its name in the dialog, when `name` is a long question. */
export interface FormBlock {name: string, fields: Field[], fix?: 'details', short?: string}
export interface FormSection {title: string, blocks: FormBlock[]}

const textOf = (c: Cell): string => (Array.isArray(c) ? c : [c]).map((r) => typeof r === 'string' ? r : r.text)
  .join('');
const lower = (s: string) => s.charAt(0).toLowerCase() + s.slice(1);
const upper = (s: string) => s.charAt(0).toUpperCase() + s.slice(1);

/** One dialog row per block, under its section: met when every field is said, else what is left to fill. */
export function formItems(sections: FormSection[]): ArriveItem[] {
  let n = 0;
  return sections.flatMap((s) => s.blocks.map((b): ArriveItem => {
    const left = b.fields.filter(([, v]) => v === null).map(([l]) => lower(l));
    const said = b.fields.filter(([, v]) => v !== null).map(([l, v]) => `${l}: ${textOf(v!)}`);
    const done = !left.length;
    // a block partly filled says what the plan filled first, then what is left to complete
    const from = b.fields.filter(([, v]) => v !== null).map(([l]) => lower(l));
    return {n: ++n, item: b.short ?? b.name, done, status: done ? 'met' : 'todo', group: s.title,
      says: done ? textOf(b.fields[0][1]!).split('\n')[0] : from.length ?
        `From the plan: ${from.join(', ')} · to complete: ${left.join(', ')}` : `To complete: ${left.join(', ')}`,
      detail: said.length ? said.join('\n') : undefined, asks: b.fields.map(([l]) => l).join(', '),
      fix: b.fix};
  }));
}

/** The document: a title, a line under it, then each section as a table of label and value (a block of several
 *  fields under its name), and the credit at the end. */
export function formBlocks(o: {title: string, sub: string, sections: FormSection[], credit: string}): Block[] {
  const out: Block[] = [{heading: o.title}, {p: o.sub, small: true}];
  for (const s of o.sections) {
    const rows: Cell[][] = [];
    for (const b of s.blocks) {
      const single = b.fields.length === 1 && b.fields[0][0] === b.name;
      if (!single)
        rows.push([{text: b.name, bold: true}, '']);
      for (const [label, value] of b.fields)
        rows.push([single ? {text: label, bold: true} : label, value ?? fill(b.short ?? label)]);
    }
    out.push({heading: s.title, level: 2}, {rows, widths: [1, 2.4]});
  }
  out.push({p: o.credit, small: true});
  return out;
}

/** What the plan says, for its documents. `protocol` the study protocol's rows as Plan shows them; `ethics` the parts
 *  of its ethics text; `perGroup` the analyzable animals of each group, the control first (null: not sized). */
export interface PlanDocFacts {
  title: string;
  meta: StudyMeta | null;
  line: string;
  protocol: {[row: string]: string};
  ethics: {design: string, animals: string, reduction: string};
  /** the study's last day, from randomization */
  end: number;
  /** the volume at randomization, mm³ */
  vRand: number;
  humane: number;
  /** the mean diameter at the humane limit, formatted (cm) */
  diameter: string;
  perGroup: {name: string, n: number | null}[];
  /** the animals to implant; null when the plan is not sized */
  implant: number | null;
  /** the measurement days */
  measurements: number;
  /** each treated arm's kind of therapy, in words */
  therapies: string[];
}

const iso = (d: Date | null) => d ? d.toISOString().slice(0, 10) : null;
const nonEmpty = (s: string | undefined) => s?.trim() ? s.trim() : null;

/** The study plan, after the blocks and labels of the ARRIVE study plan (NC3Rs). */
export function studyPlanForm(f: PlanDocFacts): FormSection[] {
  const m = f.meta;
  const a = m?.animals;
  const receives = receivesLine(m?.groups ?? []);
  // each group's size as the plan gives it, an unsized one said so: the plan sized them all, whatever it found
  const perGroup = f.perGroup.map((g) => `${g.name}: ${g.n ?? 'over 40 (no size up to 40 reaches the target power)'}`)
    .join(', ');
  const animalsMissing = !!animalsLine(a).missing.length;
  return [
    {title: 'Section 1: Study details, animals, procedures and welfare', blocks: [
      {name: 'Study details', fields: [['Study title', f.title], ['Grant code', null],
        ['Start date', m?.start ?? null], ['End date', iso(dateOf(m?.start ?? null, f.end))],
        ['Project licence or permit number', null], ['Project lead', nonEmpty(m?.runBy || m?.owner)],
        ['Protocol numbers', null],
        ['Expected severity', null], ['Primary responsible person and contact', null], ['Secondary contact', null]]},
      {name: 'Experimental animals', fix: animalsMissing ? 'details' : undefined, fields: [
        // (the atlas's models are all grown in mice)
        ['Species', nonEmpty(a?.species) ?? 'Mouse'], ['Strain or genotype', nonEmpty(a?.strain)],
        ['Sex', nonEmpty(a?.sex)],
        ['Age', nonEmpty(a?.age)], ['Weight', nonEmpty(a?.weight)], ['Source', nonEmpty(a?.source)],
        ['Number', f.protocol['Animals'] ?? null]]},
      {name: 'Experimental procedures', fix: receives ? undefined : 'details', fields: [
        ['Tumor implantation', `${f.line} tumor cells; randomization into groups when the tumors reach ${f.vRand} ` +
          'mm³ (day 0)'],
        ['Treatment by group (compound, dose, route, schedule)', receives?.detail ?? null],
        ['Measurements', f.protocol['Measurements'] ?? null], ['Duration', f.protocol['End of study'] ?? null],
        ['Surgical procedures', null], ['Anaesthesia', null], ['Analgesia', null], ['Locations', null],
        ['Acclimatisation period', null]]},
      {name: 'Animal care and monitoring', fields: [['Adverse events', null],
        ['Humane endpoints', f.protocol['Humane limit'] ?? null],
        ['Welfare monitoring', [`Tumor volume and body weight on ${f.measurements} measurement days; `,
          fill('other welfare checks')]],
        ['Clinical assessment form attached (yes or no)', null], ['Changes in housing', null],
        ['Restrictions in veterinary care', null]]},
      {name: 'Risks', fields: [['Emergency procedures', null], ['Risk to personnel', null]]},
      {name: 'Personnel', fields: [['Name, procedures, trained and competent', null]]},
    ]},
    {title: 'Section 2: Rigour', blocks: [
      {name: 'Study design and sample size', fields: [['Experimental groups', f.protocol['Design'] ?? null],
        ['Experimental unit', 'One animal'],
        ['Sample size per group', perGroup || null],
        ['Justification for sample size', f.ethics.design],
        ['Experimental Design Assistant diagram (link and access code)', null]]},
      {name: 'Inclusion and exclusion criteria', fields: [
        ['Inclusion criteria', `Animals with a tumor that has taken and is in the volume window at randomization ` +
          `(${f.vRand} mm³)`],
        ['Exclusion criteria', [`An animal whose tumor exceeds the humane limit (${f.humane} mm³) is removed; its ` +
          'measurements up to removal stay in the analysis. A measurement is left out only by a person, with the ' +
          'reason recorded. ', fill('other criteria')]],
        ['Expected attrition', f.ethics.animals]]},
      {name: 'Randomization and blinding', fields: [['Method of allocation', f.protocol['Randomization'] ?? null],
        // what the plan itself does about them; the lab's own (cage, order of measurement) left to fill
        ['Strategy to minimise confounders', ['Groups balanced on tumor volume at randomization; every group ' +
          'measured on the same days. ', fill('other confounders, e.g. cage position, order of measurement')]],
        ['Blinding strategy', m?.blinded ? ['Measurements are entered blind to the groups; the study log records ' +
          'when the groups are revealed. ', fill('who is blinded at dosing and analysis')] :
          ['Measurements are not entered blind (blinded entry is off in the study’s details). ',
            fill('who knows the groups at each step, or why not blinded')]]]},
      {name: 'Outcome measures and statistics', fields: [['Outcome measures', 'Tumor volume and body weight'],
        ['Primary outcome', 'Tumor growth rate (slope of log tumor volume) of each animal'],
        ['Analysis plans', [f.protocol['Outcome and analysis'], f.protocol['Multiple comparisons']]
          .filter((x) => x).join('. ') || null]]},
    ]},
    {title: 'Section 3: Sign-off', blocks: [
      {name: 'Sign-off', fields: [['Primary responsible person (name, date)', null],
        ['Project lead (name, date)', null], ['Internal sign-off (name, role, date)', null]]},
    ]},
  ];
}

export const STUDY_PLAN_CREDIT = 'Structure after the ARRIVE study plan, NC3Rs, arriveguidelines.org. Prepared with ' +
  'PMX Studio.';
export const NTS_CREDIT = 'Template: Commission Implementing Decision (EU) 2020/569, Annex I Part A. Prepared with ' +
  'PMX Studio.';

/** The non-technical summary's labels in the decision's own words, English and Spanish. */
const NTS = {
  en: {
    title: 'Non-technical project summary',
    sections: ['Project', 'Objectives and predicted benefits of the project', 'Predicted harms',
      'Application of the Three Rs'],
    projectTitle: 'Title of the project', duration: 'Duration of project (in months)',
    keywords: 'Key words (maximum of 5)', purpose: 'Purpose of project',
    objectives: 'Describe the objectives of the project (for example, addressing certain scientific unknowns, or ' +
      'scientific or clinical needs).',
    benefits: 'What are the potential benefits likely to derive from this project? Explain how science could be ' +
      'advanced, or humans, animals or environment may ultimately benefit from the project. Where applicable, ' +
      'differentiate between short-term benefits (within the duration of the project) and long-term benefits ' +
      '(which may accrue after the project is finished).',
    procedures: 'In what procedures will the animals typically be used (for example, injections, surgical ' +
      'procedures)? Indicate the number and duration of these procedures.',
    impacts: 'What are the expected impacts/adverse effects on the animals, for example pain, weight loss, ' +
      'inactivity/reduced mobility, stress, abnormal behavior, and the duration of those effects?',
    species: 'What species and numbers of animals are expected to be used? What are the expected severities and the ' +
      'numbers of animals in each severity category (per species)?',
    severity: 'Estimated numbers per severity (non-recovery, mild, moderate, severe)',
    fate: 'What will happen to the animals kept alive at the end of the procedure? (estimated number to be reused, ' +
      'to be returned to habitat/husbandry system, to be rehomed)',
    fateWhy: 'Please provide reasons for the planned fate of the animals after the procedure.',
    replacement: '1. Replacement. State which non-animal alternatives are available in this field and why they ' +
      'cannot be used for the purposes of the project.',
    reduction: '2. Reduction. Explain how the numbers of animals for this project were determined. Describe steps ' +
      'that have been taken to reduce the number of animals to be used, and principles used to design studies.',
    refinement: '3. Refinement. Give examples of the specific measures (e.g., increased monitoring, post-operative ' +
      'care, pain management, training of animals) to be taken, in relation to the procedures, to minimise welfare ' +
      'costs (harms) to the animals.',
    choice: 'Explain the choice of species and the related life stages.',
    retro: 'Project selected for retrospective assessment (deadline; contains severe procedures; uses non-human ' +
      'primates; other reason)',
    translational: 'Translational and applied research',
  },
  es: {
    title: 'Resumen no técnico del proyecto',
    sections: ['Proyecto', 'Objetivos y beneficios previstos del proyecto', 'Daños previstos',
      'Aplicación de la regla de las «tres erres»'],
    projectTitle: 'Título del proyecto', duration: 'Duración del proyecto (en meses)',
    keywords: 'Palabras clave (5 como máximo)', purpose: 'Finalidad del proyecto',
    objectives: 'Describa los objetivos del proyecto (por ejemplo, despejar algunas incógnitas científicas o ' +
      'atender necesidades científicas o clínicas).',
    benefits: '¿Cuáles son los probables beneficios potenciales de este proyecto? Explique la medida en que el ' +
      'proyecto puede servir al progreso científico o en que, en última instancia, puede beneficiar a los seres ' +
      'humanos, los animales o el medio ambiente. Cuando proceda, diferencie los beneficios a corto plazo (que se ' +
      'producen durante el período de ejecución del proyecto) de los beneficios a largo plazo (que pueden generarse ' +
      'una vez finalizado el proyecto).',
    procedures: '¿En qué procedimientos se utilizarán normalmente los animales (por ejemplo, inyecciones, ' +
      'procedimientos quirúrgicos)? Indique el número de procedimientos y su duración.',
    impacts: '¿Cuáles son las repercusiones o los efectos adversos previstos en los animales, por ejemplo, dolor, ' +
      'pérdida de peso, inactividad o movilidad reducida, estrés o anomalías del comportamiento, y cuál será la ' +
      'duración de esos efectos?',
    species: '¿Qué especies y qué número de animales está previsto utilizar? ¿Cuáles son los niveles de severidad ' +
      'esperados y el número de animales en cada categoría de severidad (por especies)?',
    severity: 'Cifras estimadas, por nivel de severidad (sin posibilidad de recuperación, leve, moderada, intensa)',
    fate: '¿Qué ocurrirá con los animales que se mantengan vivos al término del procedimiento? (cifra estimada de ' +
      'animales que serán reutilizados, devueltos a su hábitat/sistema de cría, realojados)',
    fateWhy: 'Le rogamos explique los motivos del destino previsto de los animales después del procedimiento.',
    replacement: '1. Reemplazo. Indique las alternativas disponibles en este campo que no requieren el uso de ' +
      'animales y el motivo por el que no pueden utilizarse para los fines del proyecto.',
    reduction: '2. Reducción. Explique cómo se ha determinado el número de animales necesario para este proyecto. ' +
      'Describa las medidas adoptadas para reducir el número de animales que van a utilizarse, así como los ' +
      'principios que inspiran los estudios.',
    refinement: '3. Refinamiento. Ofrezca ejemplos de medidas específicas (por ejemplo, intensificación del ' +
      'seguimiento, atención posoperatoria, gestión del dolor, entrenamiento de los animales) que deban adoptarse ' +
      'en relación con los procedimientos para minimizar los costes (daños) para el bienestar de los animales.',
    choice: 'Explique la elección de las especies y las etapas de vida correspondientes.',
    retro: 'Proyecto seleccionado para la evaluación retrospectiva (plazo; incluye procedimientos severos; ' +
      'utiliza primates no humanos; otros motivos)',
    translational: 'Investigación traslacional y aplicada',
  },
};
export type NtsLanguage = keyof typeof NTS;
export const ntsTitle = (lang: NtsLanguage) => NTS[lang].title;

/** The non-technical summary, Part A of the EU template: labels in `lang`, the plan's text in English. Each field is
 *  a row of the dialog, under the template's headings. */
export function ntsForm(f: PlanDocFacts, lang: NtsLanguage): FormSection[] {
  const L = NTS[lang];
  const m = f.meta;
  const receives = receivesLine(m?.groups ?? []);
  const compounds = [...new Set((m?.groups ?? []).map((g) => g.compound.trim()).filter((c) => c))];
  const treated = f.perGroup.slice(1).map((g) => g.name);
  const list = (xs: string[]) => xs.join(', ').replace(/, ([^,]*)$/, ' and $1');
  const {filled: animals, missing} = animalsLine(m?.animals);
  const species = nonEmpty(m?.animals?.species);
  const one = (short: string, label: string, value: Cell | null, fix?: 'details'): FormBlock =>
    ({name: label, short, fields: [[label, value]], fix});
  const suggest = (text: string, what: string): Run[] => [text, ' ', fill(what)];
  const title = `Efficacy of ${list(compounds.length ? compounds : treated)} in the ${f.line} tumor model`;
  const keywords = [f.line, 'tumor growth inhibition', ...new Set(f.therapies)].slice(0, 5).join(', ');
  const treatment: Run[] = receives ? [`Treatment: ${receives.detail.split('\n').join('; ')}. `] :
    [fill('treatment by group'), '. '];
  return [
    {title: L.sections[0], blocks: [
      one('Title', L.projectTitle, suggest(title, 'edit: a project may hold several studies')),
      one('Duration', L.duration, null),
      one('Key words', L.keywords, suggest(keywords, 'edit')),
      one('Purpose', L.purpose, suggest(L.translational, 'confirm')),
    ]},
    {title: L.sections[1], blocks: [one('Objectives', L.objectives, null), one('Benefits', L.benefits, null)]},
    {title: L.sections[2], blocks: [
      one('Procedures', L.procedures, [`Implantation of ${f.line} tumor cells (once). `, ...treatment,
        `Tumor and body-weight measurement on ${f.measurements} days. About ${f.end} days from randomization to ` +
        'the end of the study.'], receives ? undefined : 'details'),
      one('Adverse effects', L.impacts, [`Tumor burden up to the humane limit of ${f.humane} mm³ (a mean diameter ` +
        `of about ${f.diameter} cm), for up to ${f.end} days; body-weight loss up to ${PLAN_WEIGHT_LIMIT}%. `,
      fill('effects of the treatments and their duration')]),
      one('Species and numbers', L.species, species && f.implant !== null ? `${species}: ${f.implant} in total` :
        null, species ? undefined : 'details'),
      one('Numbers per severity', L.severity, null),
      one('Fate of the animals', L.fate, null),
      one('Reasons for their fate', L.fateWhy, null),
    ]},
    {title: L.sections[3], blocks: [
      one('Replacement', L.replacement, null),
      one('Reduction', L.reduction,
        `${f.ethics.animals} ${upper(f.ethics.reduction.replace(/^Reduction: /, ''))}`),
      one('Refinement', L.refinement, [`Humane endpoints: ${f.protocol['Humane limit'] ?? `${f.humane} mm³`}; ` +
        `body-weight loss over ${PLAN_WEIGHT_LIMIT}%. ` +
        `Tumors and body weight are measured on ${f.measurements} days. `, fill('other refinement measures')]),
      one('Choice of species', L.choice, animals ? suggest(animals, missing.length ?
        `add ${missing.join(', ')}; the reasons` : 'the reasons') : null, missing.length ? 'details' : undefined),
      one('Retrospective assessment', L.retro, null),
    ]},
  ];
}

export const ARRIVE_STUDY_PLAN_URL = 'https://arriveguidelines.org/resources/study-plan';
export const NTS_URL = 'https://eur-lex.europa.eu/eli/dec_impl/2020/569/oj';
