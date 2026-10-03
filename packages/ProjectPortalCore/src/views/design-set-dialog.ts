/** Design sets and selection sets (F-23) in the app: a design set's proposals, ranked, with the
 * ones chosen for synthesis ticked; a new design set from an open table; a new selection set from
 * the program's compounds. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {ProgramRow} from '../generated/db';
import {chooseForSynthesis, createDesignSet, createSelectionSet, DesignSetItem, designSetsFrame} from '../design-sets';
import {DESIGN_SOURCES, DesignSource, proposalColumns, proposalsFromRows} from '../model/design-sets';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const field = (label: string, control: HTMLElement, name: string): HTMLElement =>
  named(divH([span(label, 'pp-field-label'), control], 'pp-field'), name);

function select(options: {code: string, name: string}[], value = ''): HTMLSelectElement {
  const el = document.createElement('select');
  el.className = 'pp-select';
  el.append(...options.map((o) => Object.assign(document.createElement('option'), {value: o.code, textContent: o.name})));
  el.value = value || options[0]?.code || '';
  return el;
}

const text = (placeholder = ''): HTMLInputElement => Object.assign(document.createElement('input'), {className: 'pp-text', placeholder});
const failed = (e: any) => grok.shell.error(e?.message ?? String(e));

/** A design set's proposals, best rank first; ticks say chosen for synthesis. Save writes the
 * changes, with who and when. */
export function openDesignSetDialog(program: ProgramRow, item: DesignSetItem, canChoose: string[],
  onSaved?: () => void | Promise<void>): DG.Dialog {
  installStyles();
  const {set} = item;
  const t = named(document.createElement('table'), 'design-proposals');
  t.className = 'pp-card-table';
  const h = t.insertRow();
  for (const x of ['Rank', 'Label', 'Structure', set.score_name || 'Score', 'Chosen', 'Compound'])
    h.append(Object.assign(document.createElement('th'), {textContent: x}));
  const boxes = item.proposals.map((p) => {
    const tr = named(t.insertRow(), `proposal-${p.label}`);
    tr.insertCell().textContent = p.rank ? String(p.rank) : '—';
    tr.insertCell().textContent = p.label;
    tr.insertCell().textContent = p.structure;
    tr.insertCell().textContent = p.score === null || p.score === undefined ? '—' : String(Number(p.score.toPrecision(4)));
    const box = named(Object.assign(document.createElement('input'), {type: 'checkbox', checked: !!p.chosen}), `choose-${p.label}`);
    box.disabled = canChoose.length > 0;
    box.title = canChoose.join('; ') || 'Chosen for synthesis';
    tr.insertCell().append(box);
    tr.insertCell().textContent = p.compound_code ?? '';
    return {p, box};
  });
  const dlg = ui.dialog({title: `${set.name} · ${set.series_code}`});
  dlg.root.dataset.u2Name = 'design-set-dialog';
  const asTable = named(ui.button('Open as table', () => {
    grok.shell.addTableView(designSetsFrame([item], `${program.code} ${set.code}`));
    dlg.close();
  }), 'design-set-open-table');
  const scroll = divV([t]);
  scroll.style.cssText = 'max-height: 420px; overflow: auto;';
  dlg.add(divV([
    named(span(`${set.source}${set.score_name ? ` · scored by ${set.score_name}${set.higher_is_better === false ? ' (lower is better)' : ''}` : ''}` +
      `${set.description ? ` · ${set.description}` : ''}`, 'pp-muted'), 'design-set-summary'),
    scroll,
    named(divH([asTable], 'pp-actions'), 'design-set-actions'),
  ]));
  if (!canChoose.length) {
    dlg.addButton('SAVE', async () => {
      try {
        const on = boxes.filter((b) => b.box.checked && !b.p.chosen).map((b) => b.p.label);
        const off = boxes.filter((b) => !b.box.checked && b.p.chosen).map((b) => b.p.label);
        const n = (on.length ? await chooseForSynthesis(program, set.id, on, true) : 0) +
          (off.length ? await chooseForSynthesis(program, set.id, off, false) : 0);
        dlg.close();
        grok.shell.info(`${set.code}: ${n} proposal(s) changed`);
        await onSaved?.();
      } catch (e: any) {
        failed(e);
      }
    });
  }
  dlg.show({width: 760});
  if (!canChoose.length)
    named(dlg.getButton('SAVE'), 'save-choices-button');
  return dlg;
}

/** A new design set from a table open in the platform: its structure, label and score columns. */
export function openNewDesignSetDialog(program: ProgramRow, series: string[], onSaved?: () => void | Promise<void>): DG.Dialog {
  installStyles();
  const tables = grok.shell.tables;
  const table = select([{code: '', name: tables.length ? 'Choose an open table…' : 'No table is open'},
    ...tables.map((t, i) => ({code: String(i), name: `${t.name} (${t.rowCount} rows)`}))]);
  const seriesInput = select(series.map((s) => ({code: s, name: s})));
  const name = text('e.g. KIN-A para-substituents');
  const source = select(DESIGN_SOURCES.map((s) => ({code: s, name: s})), 'manual');
  const score = select([{code: '', name: 'No score'}]);
  const scoreName = text('What the score is, e.g. predicted pIC50');
  const lower = Object.assign(document.createElement('input'), {type: 'checkbox'});
  const found = named(span('', 'pp-muted'), 'design-set-columns');
  const errors = named(divV([], 'pp-warnings'), 'design-set-errors');
  const df = () => table.value ? tables[Number(table.value)] : null;
  const refresh = () => {
    const d = df();
    const cols = d ? d.columns.toList().map((c) => ({name: c.name, semType: c.semType, numeric: c.isNumerical})) : [];
    const guess = proposalColumns(cols);
    score.replaceChildren(...[{code: '', name: 'No score'}, ...cols.filter((c) => c.numeric).map((c) => ({code: c.name, name: c.name}))]
      .map((o) => Object.assign(document.createElement('option'), {value: o.code, textContent: o.name})));
    score.value = guess.score ?? '';
    found.textContent = d ? `Structures from ${guess.structure ?? '— none found —'}, labels from ${guess.label ?? 'row numbers'}` : '';
  };
  table.onchange = refresh;
  const dlg = ui.dialog({title: `New design set · ${program.code}`});
  dlg.root.dataset.u2Name = 'new-design-set-dialog';
  dlg.add(named(divV([
    field('Table', table, 'design-set-table-input'),
    found,
    field('Series', seriesInput, 'design-set-series-input'),
    field('Name', name, 'design-set-name-input'),
    field('Source', source, 'design-set-source-input'),
    field('Score', score, 'design-set-score-input'),
    field('Score is', scoreName, 'design-set-score-name-input'),
    field('Lower is better', lower, 'design-set-lower-input'),
    errors,
  ], 'pp-save-run'), 'new-design-set-form'));
  dlg.addButton('CREATE', async () => {
    try {
      const d = df();
      if (!d)
        throw new Error('Choose an open table with the proposed structures');
      const cols = proposalColumns(d.columns.toList().map((c) => ({name: c.name, semType: c.semType})), score.value || null);
      if (!cols.structure)
        throw new Error(`${d.name} has no structure column (a Molecule column, or one named smiles or structure)`);
      const rows = Array.from({length: d.rowCount}, (_, i) =>
        Object.fromEntries([cols.structure!, cols.label, cols.score].filter((c): c is string => !!c).map((c) => [c, d.get(c, i)])));
      const r = await createDesignSet(program, {name: name.value, series: seriesInput.value, source: source.value as DesignSource,
        scoreName: scoreName.value.trim() || null, higherIsBetter: !lower.checked,
        proposals: proposalsFromRows(rows, cols, `${seriesInput.value}-D`)});
      dlg.close();
      grok.shell.info(`${program.code}: design set ${r.set.code} created, ${r.proposals.length} proposals`);
      await onSaved?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'design-set-error'));
    }
  });
  dlg.show({width: 560});
  named(dlg.getButton('CREATE'), 'create-design-set-button');
  refresh();
  return dlg;
}

/** A new selection set: a name, the campaign, the rationale, the program's compounds ticked. */
export function openNewSelectionSetDialog(program: ProgramRow, campaigns: string[], compounds: string[],
  onSaved?: () => void | Promise<void>): DG.Dialog {
  installStyles();
  const name = text('e.g. K1 screen: indazole hits');
  const campaign = select([{code: '', name: 'No campaign'}, ...campaigns.map((c) => ({code: c, name: c}))]);
  const rationale = Object.assign(document.createElement('textarea'), {rows: 3, placeholder: 'Why these compounds'});
  const boxes = compounds.map((c) => ({c, box: named(Object.assign(document.createElement('input'), {type: 'checkbox'}), `select-${c}`)}));
  const errors = named(divV([], 'pp-warnings'), 'selection-set-errors');
  const list = divV(boxes.map((b) => divH([b.box, span(b.c)], 'pp-unit')));
  list.style.cssText = 'max-height: 240px; overflow: auto;';
  const dlg = ui.dialog({title: `New selection set · ${program.code}`});
  dlg.root.dataset.u2Name = 'new-selection-set-dialog';
  dlg.add(named(divV([
    field('Name', name, 'selection-set-name-input'),
    field('From campaign', campaign, 'selection-set-campaign-input'),
    field('Rationale', rationale, 'selection-set-rationale-input'),
    span('Compounds', 'pp-heading'),
    named(list, 'selection-set-compounds'),
    errors,
  ], 'pp-save-run'), 'new-selection-set-form'));
  dlg.addButton('CREATE', async () => {
    try {
      const r = await createSelectionSet(program, {name: name.value, campaign: campaign.value || null, rationale: rationale.value,
        compounds: boxes.filter((b) => b.box.checked).map((b) => b.c)});
      dlg.close();
      grok.shell.info(`${program.code}: selection set ${r.code} created`);
      await onSaved?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'selection-set-error'));
    }
  });
  dlg.show({width: 520});
  named(dlg.getButton('CREATE'), 'create-selection-set-button');
  return dlg;
}
