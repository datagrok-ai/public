/** Design sets and selection sets (F-23) in the app: a design set's proposals, ranked, with the
 * ones chosen for synthesis ticked; a new design set from an open table; a new selection set from
 * the program's compounds. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import {button, divH, divV, span} from '@datagrok-libraries/u2';
import {ProgramRow} from '../generated/db';
import {chooseForSynthesis, createDesignSet, createSelectionSet, DesignSetItem, designSetsFrame} from '../design-sets';
import {DESIGN_SOURCES, DesignSource, proposalColumns, proposalsFromRows} from '../model/design-sets';
import {area, check, choice, choiceValue, ErrorRegion, named, openDialog, PortalDialog, table,
  text} from './u2-controls';

/** A design set's proposals, best rank first; ticks say chosen for synthesis. Save writes the
 * changes, with who and when. */
export function openDesignSetDialog(program: ProgramRow, item: DesignSetItem, canChoose: string[],
  onSaved?: () => void | Promise<void>): PortalDialog {
  const {set} = item;
  let changed = 0;
  const title = `${set.name} · ${set.series_code}`;
  const dlg: PortalDialog = openDialog({title, name: 'design-set-dialog', width: 760},
    () => {
      const boxes = item.proposals.map((p) => ({p, box: check(`choose-${p.label}`, undefined, !!p.chosen,
        {enabled: canChoose.length === 0, tooltipText: canChoose.join('; ') || 'Chosen for synthesis'})}));
      const t = table('design-proposals', [
        {header: 'Rank', render: (b: typeof boxes[number]) => b.p.rank ? String(b.p.rank) : '—'},
        {header: 'Label', render: (b) => b.p.label},
        {header: 'Structure', render: (b) => b.p.structure},
        {header: set.score_name || 'Score', render: (b) =>
          b.p.score === null || b.p.score === undefined ? '—' : String(Number(b.p.score.toPrecision(4)))},
        {header: 'Chosen', render: (b) => b.box.root},
        {header: 'Compound', render: (b) => b.p.compound_code ?? ''},
      ], boxes, {name: (b) => `proposal-${b.p.label}`});
      const asTable = named(button('Open as table', () => {
        grok.shell.addTableView(designSetsFrame([item], `${program.code} ${set.code}`));
        dlg.close();
      }), 'design-set-open-table');
      const scroll = divV([t.root]);
      scroll.style.cssText = 'max-height: 420px; overflow: auto;';
      return {
        content: divV([
          named(span(`${set.source}${set.score_name ? ` · scored by ${set.score_name}` +
            `${set.higher_is_better === false ? ' (lower is better)' : ''}` : ''}` +
            `${set.description ? ` · ${set.description}` : ''}`, 'pp-muted'), 'design-set-summary'),
          scroll,
          named(divH([asTable], 'pp-actions'), 'design-set-actions'),
        ]),
        buttons: canChoose.length ? [] : [{label: 'SAVE', name: 'save-choices-button', primary: true,
          run: async () => {
            const on = boxes.filter((b) => b.box.value.value && !b.p.chosen).map((b) => b.p.label);
            const off = boxes.filter((b) => !b.box.value.value && b.p.chosen).map((b) => b.p.label);
            changed = (on.length ? await chooseForSynthesis(program, set.id, on, true) : 0) +
              (off.length ? await chooseForSynthesis(program, set.id, off, false) : 0);
          },
          then: async () => {
            grok.shell.info(`${set.code}: ${changed} proposal(s) changed`);
            await onSaved?.();
          }}],
      };
    });
  return dlg;
}

/** A new design set from a table open in the platform: its structure, label and score columns. */
export function openNewDesignSetDialog(program: ProgramRow, series: string[],
  onSaved?: () => void | Promise<void>): PortalDialog {
  const tables = grok.shell.tables;
  const errors = new ErrorRegion('design-set-errors', 'design-set-error');
  let created = '';
  const title = `New design set · ${program.code}`;
  return openDialog({title, name: 'new-design-set-dialog', width: 560, errors}, () => {
    const found = named(span('', 'pp-muted'), 'design-set-columns');
    const table = choice('design-set-table-input', 'Table',
      tables.map((t, i) => ({code: String(i), name: `${t.name} (${t.rowCount} rows)`})), '',
      {blank: tables.length ? 'Choose an open table…' : 'No table is open', onChanged: () => refresh()});
    const seriesInput = choice('design-set-series-input', 'Series', series.map((s) => ({code: s, name: s})));
    const name = text('design-set-name-input', 'Name', {placeholder: 'e.g. KIN-A para-substituents'});
    const sources = DESIGN_SOURCES.map((s) => ({code: s, name: s}));
    const source = choice('design-set-source-input', 'Source', sources, 'manual');
    const score = choice('design-set-score-input', 'Score', [], '', {blank: 'No score'});
    const scoreName = text('design-set-score-name-input', 'Score is',
      {placeholder: 'What the score is, e.g. predicted pIC50'});
    const lower = check('design-set-lower-input', 'Lower is better');
    const df = () => choiceValue(table) ? tables[Number(choiceValue(table))] : null;
    const refresh = () => {
      const d = df();
      const cols = d ? d.columns.toList().map((c) => ({name: c.name, semType: c.semType, numeric: c.isNumerical})) : [];
      const guess = proposalColumns(cols);
      const numeric = cols.filter((c) => c.numeric).map((c) => ({value: c.name, label: c.name}));
      score.setItems([{value: '', label: 'No score'}, ...numeric]);
      score.value.value = guess.score ?? '';
      found.textContent = d ? `Structures from ${guess.structure ?? '— none found —'}, ` +
        `labels from ${guess.label ?? 'row numbers'}` : '';
    };
    refresh();
    return {
      content: named(divV([table.root, found, seriesInput.root, name.root, source.root, score.root, scoreName.root,
        lower.root, errors.root], 'pp-save-run'), 'new-design-set-form'),
      buttons: [{label: 'CREATE', name: 'create-design-set-button', primary: true,
        run: async () => {
          const d = df();
          if (!d)
            throw new Error('Choose an open table with the proposed structures');
          const cols = proposalColumns(d.columns.toList().map((c) => ({name: c.name, semType: c.semType})),
            choiceValue(score) || null);
          if (!cols.structure)
            throw new Error(`${d.name} has no structure column (a Molecule column, or one named smiles or structure)`);
          const keep = [cols.structure!, cols.label, cols.score].filter((c): c is string => !!c);
          const rows = Array.from({length: d.rowCount},
            (_, i) => Object.fromEntries(keep.map((c) => [c, d.get(c, i)])));
          const r = await createDesignSet(program, {name: name.value.value, series: choiceValue(seriesInput),
            source: choiceValue(source) as DesignSource, scoreName: scoreName.value.value.trim() || null,
            higherIsBetter: !lower.value.value,
            proposals: proposalsFromRows(rows, cols, `${choiceValue(seriesInput)}-D`)});
          created = `${program.code}: design set ${r.set.code} created, ${r.proposals.length} proposals`;
        },
        then: async () => {
          grok.shell.info(created);
          await onSaved?.();
        }}],
    };
  });
}

/** A new selection set: a name, the campaign, the rationale, the program's compounds ticked. */
export function openNewSelectionSetDialog(program: ProgramRow, campaigns: string[], compounds: string[],
  onSaved?: () => void | Promise<void>): PortalDialog {
  const errors = new ErrorRegion('selection-set-errors', 'selection-set-error');
  let created = '';
  const title = `New selection set · ${program.code}`;
  return openDialog({title, name: 'new-selection-set-dialog', width: 520, errors},
    () => {
      const name = text('selection-set-name-input', 'Name', {placeholder: 'e.g. K1 screen: indazole hits'});
      const campaignOptions = campaigns.map((c) => ({code: c, name: c}));
      const campaign = choice('selection-set-campaign-input', 'From campaign', campaignOptions, '',
        {blank: 'No campaign'});
      const rationale = area('selection-set-rationale-input', 'Rationale',
        {placeholder: 'Why these compounds', rows: 3});
      const boxes = compounds.map((c) => ({c, box: check(`select-${c}`, c)}));
      const list = divV(boxes.map((b) => b.box.root), 'pp-checklist');
      list.style.cssText = 'max-height: 240px; overflow: auto;';
      return {
        content: named(divV([name.root, campaign.root, rationale.root, span('Compounds', 'pp-heading'),
          named(list, 'selection-set-compounds'), errors.root], 'pp-save-run'), 'new-selection-set-form'),
        buttons: [{label: 'CREATE', name: 'create-selection-set-button', primary: true,
          run: async () => {
            const r = await createSelectionSet(program, {name: name.value.value,
              campaign: choiceValue(campaign) || null,
              rationale: rationale.value.value, compounds: boxes.filter((b) => b.box.value.value).map((b) => b.c)});
            created = `${program.code}: selection set ${r.code} created`;
          },
          then: async () => {
            grok.shell.info(created);
            await onSaved?.();
          }}],
      };
    });
}
