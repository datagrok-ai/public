/* Compare: studies with measurements side by side, each treated group's effect against the patients' exposure.
   What each study kept is read as saved (its groups, its NCA and patients' exposure, and the analysis's TGI on its
   evaluation day, pmx.exposure); a study analyzed before that was kept gives 1 − T/C of the median volumes against its
   own control (exposure.ts medianTgi), said so. It describes what was measured: no direction, no cause, nothing
   beyond the groups. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {HumanExposure, HumanMultiple, cvOk, drugMultiples, fuOk, humanMultiple, humansInto, medianTgi, multipleText,
  parseExposure, patientsRange, savedEffect, sourceText} from './exposure';
import type {ArmExposure, CurveRow} from './exposure';
import {PreviousStudy, controlOf, readDayOf} from './study-effect';
import {armParts, detailsName, metaOf} from './study-meta';
import {describeArm} from './arm-name';
import {previousOfSaved} from './previous-studies';
import {entryText, open, save, saveError} from './store';
import {EXPOSURE_TAG, exposureOfTag, parseStudy} from './study-file';
import {humanExposureDialog} from './human-exposure';
import type {LibraryEntry} from './library';
import {boundText, exposureSection, sectionElement} from './translation-report';
import {reportHtml} from './report';
import {esc, printSheet} from './sheet';
import {initQuietly} from './band-chart';
import {STANDS_ALONE} from './stages';

const NO_MOUSE = 'the study kept no mouse exposure (NCA)';
/** Under a combination's multiples, and in the report: whose AUC they use. */
export const COMBO_AUC = 'each drug at its AUC in this combination group\'s mice';

/** One treated group of one study, as Compare shows it. */
export interface CompareRow {
  study: string, model: string, group: string,
  /** the drug (a combination's drugs joined by ' + ') and its key ('' for a combination: never placed) */
  drug: string, key: string, treatment: string, dose: string | null,
  /** TGI (%), on `day`; `analyzed`: the analysis's (1 − ΔT/ΔC), else 1 − T/C of the medians */
  tgi: number | null, day: number | null, analyzed: boolean,
  multiple: HumanMultiple | null,
  /** where the patients' exposure came from; null: not given */
  source: string | null, cv: number | null, unbound: boolean,
  /** the group's NCA exposure and its drug's patients' exposure, as kept (for the report) */
  x?: ArmExposure, human?: HumanExposure,
  /** a combination with NCA's exposure kept for its drugs: each drug against its patients' exposure (not charted) */
  drugs?: {drug: string, m: HumanMultiple | null, x?: ArmExposure, human?: HumanExposure}[],
}

/** A study's groups, in the order of its rows. */
const groupsIn = (s: PreviousStudy) => [...new Set(s.rows.map((r) => r.group).filter((g): g is string => !!g))];

/** A study's treated groups: each with the analysis's TGI when the study kept it (else 1 − T/C of the medians on
 *  the last day with 3 animals in it and its control), and its exposure as a multiple of the patients'. */
export function compareRows(s: PreviousStudy, model: string): CompareRow[] {
  const groups = groupsIn(s);
  const control = controlOf(groups, s.meta);
  const rows = s.rows.filter((r): r is CurveRow => r.group !== null && r.day !== null && r.volume !== null);
  const {humans, nca} = parseExposure(s.exposure);
  const effect = savedEffect(s.exposure);
  return groups.filter((g) => g !== control).map((g) => {
    const parts = armParts(g, s.meta);
    const one = parts.length === 1 ? parts[0] : null;
    const h = one ? humans.get(one.key) : undefined;
    const x = nca?.arms.get(g)?.exposure;
    const multiple = one && h && h.auc !== null ? nca && (x || nca.pk.get(g)) ?
      humanMultiple(g, x, nca.pk.get(g) ?? null, humans, s.meta) : {why: NO_MOUSE} : null;
    const analyzed = effect?.tgi.get(g);
    // a combination's source: each of its drugs that has one, named (its multiples are each drug's)
    const given = parts.flatMap((p) => humans.get(p.key) ? [`${p.drug}: ${sourceText(humans.get(p.key)!)}`] : []);
    return {study: s.name, model, group: g, drug: parts.map((p) => p.drug).join(' + ') || g, key: one?.key ?? '',
      treatment: describeArm(parts) || g, dose: one?.dose ? `${one.dose.value} ${one.dose.unit}` : null,
      ...analyzed !== undefined ? {tgi: analyzed, day: effect!.day, analyzed: true} :
        {tgi: control ? medianTgi(rows, g, control, null) : null, day: control ? readDayOf(s, g, control) : null,
          analyzed: false},
      multiple, source: one ? h ? sourceText(h) : null : given.length ? given.join('; ') : null,
      cv: h && cvOk(h.cv) ? h.cv : null,
      unbound: !!h && fuOk(h.fuMouse) && fuOk(h.fuHuman), x: nca?.arms.get(g), human: h,
      ...parts.length > 1 && nca?.drugs?.get(g) ?
        {drugs: drugMultiples(g, nca.drugs.get(g)!, humans, s.meta).map(({drug, key, m}) => ({drug, m,
          x: nca.drugs!.get(g)!.find((d) => d.key === key)?.x, human: humans.get(key)}))} : {}};
  });
}

/** The rows the chart shows: those with a TGI and a multiple of the patients' exposure. */
export const chartRows = (rows: CompareRow[]) =>
  rows.filter((r) => r.tgi !== null && Number.isFinite(r.tgi) && r.multiple?.ratio !== undefined);

/** The one line when the studies are of more than one model; null: one model. */
export const modelWarning = (models: string[]) => new Set(models.filter((m) => m)).size > 1 ?
  'Different models: tumour growth inhibition is not directly comparable across models.' : null;

/** One hint per drug no study gives the patients' exposure of (its first study, where it would be kept). */
export function hintRows(rows: CompareRow[]): {drug: string, key: string, study: string}[] {
  const given = new Set(rows.filter((r) => r.source).map((r) => r.key));
  const out = new Map<string, {drug: string, key: string, study: string}>();
  for (const r of rows) {
    if (r.key && !given.has(r.key) && !out.has(r.key))
      out.set(r.key, {drug: r.drug, key: r.key, study: r.study});
  }
  return [...out.values()];
}

/** What a TGI is, when it is not the analysis's (shown with an asterisk). */
export const tgiNote = (r: CompareRow) => r.analyzed ? null :
  `1 − T/C of medians${r.day !== null ? `, day ${r.day}` : ''} (not from the analysis)`;

/** The footnote of the TGIs marked with an asterisk; null: every TGI shown is the analysis's. */
export const tgiFootnote = (rows: CompareRow[]): string | null => rows.some((r) => !r.analyzed && r.tgi !== null) ?
  '* The study kept no analysis: 1 − T/C of the median volumes, on the last day with 3 animals in the group and ' +
  'its control.' : null;

const pct = (x: number) => `${Math.round(x)}%`;
/** A multiple as the table shows it: its arrow, and averaged over the week when no dose is given at weekends (the
 *  daily AUC is then the dose's times 5/7, as Analyze says: the AUC per dose alone would not give it). */
export const shownMultiple = (m: HumanMultiple & {ratio: number}) => multipleText(m.ratio) +
  (m.range === 'above' ? ' ↑' : m.range === 'below' ? ' ↓' : '') + (m.weekly ? ', averaged over the week' : '');
const rangeWord = (m: HumanMultiple) => m.range === 'above' ? ', above the patients\' range' :
  m.range === 'below' ? ', below the patients\' range' : '';

/** One clause per drug on what its groups gave, from the numbers alone (the study named when the drug is in more
 *  than one); then the drugs not placed, and why. null: nothing to say. */
export function compareSentence(rows: CompareRow[]): string | null {
  const placed = chartRows(rows);
  const drugs = [...new Set(placed.map((r) => r.drug))];
  const said = drugs.map((d) => {
    const gs = placed.filter((r) => r.drug === d).sort((a, b) => a.multiple!.ratio! - b.multiple!.ratio!);
    const studies = new Set(rows.filter((r) => r.drug === d).map((r) => r.study)).size > 1;
    const each = gs.map((r) => `${pct(r.tgi!)}${r.dose ? ` at ${r.dose}` : ''}${studies ? ` in ${r.study}` : ''} ` +
      `(${multipleText(r.multiple!.ratio!)}${rangeWord(r.multiple!)}` +
      `${r.multiple!.weekly ? ', averaged over the week' : ''})`);
    const list = each.length > 1 ? `${each.slice(0, -1).join(', ')} and ${each[each.length - 1]}` : each[0];
    return `${d} gave a TGI of ${list}`;
  });
  const placedKeys = new Set(placed.map((r) => r.key));
  const noMouse = [...new Set(rows.filter((r) => r.key && !placedKeys.has(r.key) && r.multiple?.why === NO_MOUSE)
    .map((r) => r.drug))].map((d) => `${d} has no mouse exposure`);
  const missing = hintRows(rows).map((h) => `${h.drug} has no patients' exposure yet`);
  const out = [...said, ...noMouse, ...missing];
  return out.length ? `${out.join('; ')}.`.replace(/^./, (c) => c.toUpperCase()) : null;
}

/** What the exposures are and where they come from. */
export function compareFooter(rows: CompareRow[]): string {
  // (a combination's source, its drugs' each named, stays on its row: here it would repeat the single drugs')
  const sources = [...new Set(rows.filter((r) => r.source && r.key).map((r) => r.source!))];
  return `${boundText(rows.map((r) => r.human))} Mouse AUC from each study's NCA` +
    `${sources.length ? `; patients' exposure from ${sources.join('; ')}` : ''}.`;
}

/** A study file's text with the patients' exposure set to `humans`, all else as the file has it now (read again
 *  just before writing: the study may have been saved since Compare read it). */
export function withHumans(text: string, humans: Map<string, HumanExposure>): string {
  const raw = JSON.parse(text);
  raw.exposure = humansInto(raw.exposure, humans);
  return JSON.stringify(raw);
}

/** The study's tables open now (its ID and its name: a copy of the table is left alone) given the same patients'
 *  exposure as its file: a Save from Monitor then keeps it, and Analyze takes it in (analyze-view keepExposure). */
export function humansToOpen(id: string, name: string, humans: Map<string, HumanExposure>,
  tables: DG.DataFrame[] = grok.shell.tables): void {
  for (const df of id ? tables.filter((t) => metaOf(t)?.id === id && t.name === name) : [])
    df.setTag(EXPOSURE_TAG, JSON.stringify(humansInto(exposureOfTag(df), humans)));
}

/** The text of each entry: the library lists studies without reading them, so each is read now from its file; one
 *  that cannot be read (deleted, damaged) is skipped, by name. */
export async function compareTexts(entries: LibraryEntry[],
  read: (f: DG.FileInfo) => Promise<string> = (f) => grok.dapi.files.readAsText(f)):
Promise<{texts: Map<LibraryEntry, string>, skipped: string[]}> {
  const read1 = async (e: LibraryEntry): Promise<string | null> => {
    try {
      const text = e.text || (e.file ? await read(e.file) : e.row ? await entryText(e) : '');
      parseStudy(text);
      return text;
    } catch {
      return null;
    }
  };
  const texts = await Promise.all(entries.map(read1));
  return {texts: new Map(entries.flatMap((e, i): [LibraryEntry, string][] =>
    texts[i] === null ? [] : [[e, texts[i]!]])),
  skipped: entries.filter((_, i) => texts[i] === null).map((e) => e.name || e.id || 'A study')};
}

/** Each study's name as Compare shows it: with its ID when another study shown has the same name, and "running, day
 *  n" when it is still followed (its results are not final); unique, so the chart keys its series by it. */
export function studyLabels(entries: Pick<LibraryEntry, 'name' | 'id' | 'stage' | 'day'>[]): string[] {
  const count = new Map<string, number>();
  for (const e of entries)
    count.set(e.name, (count.get(e.name) ?? 0) + 1);
  const out = entries.map((e) => {
    const more = [...(count.get(e.name) ?? 0) > 1 && e.id ? [e.id] : [],
      ...e.stage === 2 ? [`running, day ${e.day ?? 0}`] : []];
    return `${e.name || e.id || 'A study'}${more.length ? ` (${more.join(', ')})` : ''}`;
  });
  // (the same name and no ID to tell them apart: numbered)
  return out.map((l, i) => out.indexOf(l) === i ? l : `${l} #${out.slice(0, i + 1).filter((x) => x === l).length}`);
}

/** The report's section for the compared studies: the sentence, the groups, the models' warning, what it does
 *  not say and the sources (translation-report.ts). */
export const compareSection = (rows: CompareRow[], models: string[]) => exposureSection({sentence:
  compareSentence(rows), warning: modelWarning(models), curve: false, rows: rows.flatMap((r) => {
    const common = {tgi: r.tgi, study: r.study, model: r.model, ...r.analyzed ? {} : {medians: {day: r.day}}};
    // a combination: a row per drug, "A + B (A)", with the multiple the table shows (as Analyze's report does)
    return r.drugs?.some((d) => d.m) ? r.drugs.map((d) => ({...common, group: `${r.treatment} (${d.drug})`, x: d.x,
      multiple: d.m, human: d.human})) : [{...common, group: r.treatment, x: r.x, multiple: r.multiple,
      human: r.human}];
  })});

/** The chart: each group at its TGI against its multiple of the patients' exposure, on a log axis ticked at each
 *  doubling, one colour per study (its name in the legend, above the plot), each drug's patients' range shaded
 *  behind. echarts, as the band chart (the platform's scatter plot gave no log ticks at these few points). */
function compareChart(shown: CompareRow[], y: string): {root: HTMLElement, dispose: () => void} {
  const root = ui.div([], {style: {height: '320px', width: '100%'}});
  const chart = initQuietly(root);
  const studies = [...new Set(shown.map((r) => r.study))];
  const ranges = [...new Set(shown.map((r) => r.cv))].map((cv) => ({cv, ...patientsRange(cv),
    drugs: [...new Set(shown.filter((r) => r.cv === cv).map((r) => r.drug))].join(', ')}));
  const xs = [...shown.map((r) => r.multiple!.ratio!), ...ranges.flatMap((b) => [b.lo, b.hi])];
  const fmt = (v: number) => `${Number(v.toPrecision(2))}×`;
  const [lo, hi] = [Math.floor(Math.log2(Math.min(...xs))), Math.ceil(Math.log2(Math.max(...xs)))];
  // a point in the right third of the axis is labelled on its left: on its right the label ran off the chart
  const leftOf = (ratio: number) => (Math.log2(ratio) - lo) / Math.max(hi - lo, 1) > 2 / 3;
  chart.setOption({
    animation: false,
    grid: {left: 56, right: 32, top: 40, bottom: 44},
    legend: {top: 0, data: studies},
    // (echarts reads the tooltip as HTML: the study's and the group's names are escaped)
    tooltip: {trigger: 'item', formatter: (p: any) => `${esc(p.seriesName)}<br>${esc(p.name)}: ` +
      `TGI ${Math.round(p.value[1])}% at ${fmt(p.value[0])} the patients' exposure`},
    xAxis: {type: 'log', logBase: 2, min: 2 ** lo, max: 2 ** hi, name: 'Multiple of the patients\' exposure',
      nameLocation: 'middle', nameGap: 28, axisLabel: {formatter: fmt}},
    yAxis: {type: 'value', name: y, nameLocation: 'middle', nameGap: 40},
    series: studies.map((s, i) => ({type: 'scatter', name: s, symbolSize: 10,
      color: DG.Color.toHtml(DG.Color.getCategoricalColor(i)),
      data: shown.filter((r) => r.study === s).map((r) => ({name: `${r.drug} ${r.dose ?? ''}`.trim(),
        value: [r.multiple!.ratio!, r.tgi!], ...leftOf(r.multiple!.ratio!) ? {label: {position: 'left'}} : {}})),
      label: {show: true, position: 'right', formatter: '{b}', fontSize: 10},
      // two groups at nearly one point (one dose in two studies): their labels stacked, not drawn over each other
      labelLayout: {moveOverlap: 'shiftY'},
      // the patients' ranges, behind the points (on the first series: drawn once)
      ...i === 0 ? {markArea: {silent: true, itemStyle: {color: 'rgba(154, 165, 177, 0.18)'},
        label: {position: 'insideTop', color: '#717581', fontSize: 10},
        data: ranges.map((b) => [{xAxis: b.lo, name: `${b.drugs}: patients' range` +
          (cvOk(b.cv) ? ` (CV ${b.cv}%)` : ' (3-fold convention)')}, {xAxis: b.hi}])}} : {}})),
  });
  const resized = ui.onSizeChanged(root).subscribe(() => chart.resize());
  return {root, dispose: () => {
    resized.unsubscribe();
    chart.dispose();
  }};
}

/** The Compare view of the selected studies with measurements (library entries, each read from its file: the
 *  library lists them without their text). */
export async function compareView(entries: LibraryEntry[],
  read?: (f: DG.FileInfo) => Promise<string>): Promise<DG.View> {
  const {texts, skipped} = await compareTexts(entries, read);
  const shownEntries = entries.filter((e) => texts.has(e));
  const view = DG.View.create();
  view.name = 'PMX Studio · Compare';
  // no Toolbox or Browse beside it (empty for this page): it stands alone, as the Overview does
  view.root.dataset[STANDS_ALONE] = 'true';
  const host = ViewHost.of(view);
  const muted = (t: string) => {
    const el = ui.divText(t);
    el.style.color = 'var(--grey-4)';
    return el;
  };
  // the chart drawn last, let go when it is drawn again or the view closes
  let disposeChart = () => {};
  host.own(() => disposeChart());
  const labels = studyLabels(shownEntries);
  const draw = () => {
    const studies = shownEntries.map((e, i) => ({e, s: previousOfSaved(parseStudy(texts.get(e)!)), label: labels[i]}));
    // each study by its label (unique: two studies of one name told apart by their IDs), in the chart's series, the
    // table, the sentence and the report
    const rows = studies.flatMap(({e, s, label}) => compareRows(s, e.line).map((r) => ({...r, study: label})));
    const sentence = compareSentence(rows);
    const warning = modelWarning(shownEntries.map((e) => e.line));
    // the chart: TGI against the multiple, log X, the patients' range shaded (each range with its drugs)
    const shown = chartRows(rows);
    disposeChart();
    let chart: HTMLElement | null = null;
    if (shown.length) {
      const y = shown.every((r) => r.analyzed) ? 'TGI from the analysis (%)' : shown.some((r) => r.analyzed) ?
        'TGI from the analysis, * 1 − T/C of medians (%)' : 'TGI, * 1 − T/C of medians (%)';
      const drawn = compareChart(shown, y);
      chart = drawn.root;
      disposeChart = drawn.dispose;
    }
    // built under the view's host: draw() also runs after an await (a patients' exposure added), when nothing builds
    // the page's width shared by the columns, long text wrapped in its cell: a long source would widen the table past
    // the window
    const wrap = (t: string) => ui.div([t], {style: {whiteSpace: 'normal', overflowWrap: 'anywhere'}});
    const table = host.build(() => new BasicTable<CompareRow>({items: rows, fixedLayout: true,
      keyOf: (r) => `${r.study}\u0000${r.group}`, columns: [
        {header: 'Study', render: (r) => wrap(r.study)},
        {header: 'Model', render: (r) => r.model},
        {header: 'Treatment', render: (r) => wrap(r.treatment)},
        {header: 'TGI from the analysis', render: (r) => {
          if (r.tgi === null)
            return '—';
          const el = ui.span([`${pct(r.tgi)}${r.analyzed ? '' : ' *'}${r.day !== null ? ` · day ${r.day}` : ''}`]);
          const note = tgiNote(r);
          if (note)
            ui.tooltip.bind(el, note);
          return el;
        }},
        {header: '× patients', render: (r) => {
          // a combination: each drug's multiple, its name first
          if (r.drugs?.some((d) => d.m)) {
            // each drug at the AUC NCA read in this group's own mice: not the same as the drug's group alone (1.9×
            // here, 2.2× there is two groups' exposure, not two computations)
            return ui.divV([...r.drugs.map((d) => {
              const el = ui.divText(`${d.drug} ${d.m?.ratio !== undefined ? shownMultiple(d.m) : '—'}`);
              if (d.m)
                ui.tooltip.bind(el, d.m.ratio !== undefined ? d.m.tip : d.m.why);
              return el;
            }), muted(COMBO_AUC)]);
          }
          const m = r.multiple;
          if (m?.ratio !== undefined) {
            const el = ui.span([shownMultiple(m)]);
            ui.tooltip.bind(el, m.tip);
            return el;
          }
          return ui.divV([ui.span(['—']), muted(m?.why ?? (r.key ? 'no patients\' exposure' : 'a combination'))]);
        }},
        {header: 'Patients\' exposure from', render: (r) => wrap(r.source ? `${r.source}${r.cv !== null ?
          ` · CV ${r.cv}%` : ''}` : 'Not given')},
      ]}));
    // one hint per drug without the patients' exposure: typed here, kept with its study
    const hints = hintRows(rows).map((hint) => {
      const add = ui.button('Add patients\' exposure…', () => {
        const {e, s} = studies.find((x) => x.label === hint.study)!;
        const {file, row} = e;
        if (!file && !row)
          return grok.shell.warning('This study is not saved in your files: open it in Analyze to keep it');
        const names = groupsIn(s).map((g) => detailsName(g, s.meta));
        // the study as it is now, changed and written back: its file, or its row at the version just read
        const update = async (change: (now: string) => string): Promise<string> => {
          if (file) {
            const text = change(await grok.dapi.files.readAsText(file));
            await grok.dapi.files.writeAsText(file, text);
            return text;
          }
          const {doc, saved} = await open('study', row!.id);
          const text = change(JSON.stringify(doc));
          await save('study', parseStudy(text), saved);
          return text;
        };
        void humanExposureDialog(names, parseExposure(s.exposure).humans, (m) => {
          void update((now) => withHumans(now, m)).then((text) => {
            texts.set(e, text);
            humansToOpen(e.id, e.name, m);
            draw();
          }).catch((err) => grok.shell.error(row ? saveError(err, 'study') :
            `The patients' exposure could not be saved with ${esc(hint.study)}: ${esc(err)}`));
        }).catch((err) => grok.shell.error(String(err)));
      });
      return ui.divH([muted(`${hint.drug}: add the patients' exposure (a projection is fine) to place it on the ` +
        'chart.'), add], {style: {alignItems: 'center', gap: '8px'}});
    });
    // the report: the chart as drawn, then the section (what is on screen besides the chart is in it)
    const report = ui.button('Report', () => {
      const section = compareSection(rows, shownEntries.map((e) => e.line));
      void reportHtml({title: 'Compare', stage: 'Compare', meta: null, content: page, annex: section ?
        [['Exposure and patients', sectionElement(section)]] : []}).then(printSheet,
        (err) => grok.shell.error(`The report could not be built: ${err}`));
    });
    const noprint = (els: HTMLElement[]) => els.map((el) => {
      el.classList.add('pmx-noprint');
      return el;
    });
    const footnote = tgiFootnote(rows);
    const page = ui.divV([...noprint([report, ...sentence ? [ui.divText(sentence)] : [],
      ...warning ? [muted(warning)] : [], ...skipped.map((n) => muted(`${n} could not be read (moved, deleted or ` +
        'damaged): left out.'))]), ...chart ? [chart] : [], ...noprint([table.root,
      ...footnote ? [muted(footnote)] : [], ...hints, muted(compareFooter(rows))])],
    {style: {gap: '12px', padding: '16px', overflow: 'auto'}});
    view.root.replaceChildren(page);
  };
  draw();
  return view;
}
