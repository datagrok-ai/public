import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {badge, button, Control, divH, divV, icon, Section, span} from '@datagrok-libraries/u2';
import {appView, viewers} from '@datagrok-libraries/u2/src/dg/index.js';

import {JmpTable} from './jmp-types';
import {addSavedAnalyses, dataColumns, jmpTableToDataFrame, savedAnalyses, skippedColumns} from './jmp-to-dataframe';
import {JmpScriptConversion} from './jsl-to-viewers';

/** File-browser preview of a JMP table (JMP_PLUGIN.md §6.1): what the file is, the columns that
 * won't come across, the saved analyses and what they open as, the data, and Open in the ribbon. */
export function createJmpPreview(table: JmpTable, fileName: string): DG.ViewBase {
  const tableName = fileName.replace(/\.jmp$/i, '');
  const df = jmpTableToDataFrame(table, tableName);
  const skipped = skippedColumns(table);
  const analyses = savedAnalyses(table);

  const content = Control.build(() => {
    const summary = divH([
      badge(`JMP ${table.info.jmpVersion}`, {variant: 'accent'}),
      span(`${table.info.nrows.toLocaleString()} rows × ${dataColumns(table).length} columns`, 'jmp-preview-size'),
      span(`Saved ${table.info.savetime.toLocaleString()}`, 'jmp-preview-saved'),
    ], 'jmp-preview-summary');
    summary.dataset.u2Name = 'summary';
    Object.assign(summary.style, {gap: '12px', alignItems: 'center', padding: '6px 8px', flex: 'none'});

    const parts: (HTMLElement | Control)[] = [summary];
    if (skipped.length > 0) {
      const section = new Section({
        title: 'Columns not imported',
        hint: `${skipped.length} of ${table.info.ncols}`,
        expanded: false,
        persistKey: 'jmp.preview.skipped',
      });
      section.name = 'skipped columns';
      section.add(divV(skipped.map((s) => {
        const row = divH([span(s.name, 'jmp-skipped-name'), span(s.reason, 'jmp-skipped-reason')]);
        row.style.gap = '12px';
        row.dataset.u2Item = s.name;
        return row;
      })));
      section.root.style.flex = 'none';
      parts.push(section);
    }

    if (analyses.length > 0)
      parts.push(analysesSection(analyses));

    const grid = viewers.grid(df);
    grid.name = 'data';
    Object.assign(grid.root.style, {flex: '1 1 0', minHeight: '0', width: '100%'});
    parts.push(grid);

    const root = divV(parts, 'jmp-preview');
    Object.assign(root.style, {width: '100%', height: '100%'});
    return root;
  });

  const converted = analyses.filter((a) => a.viewers.length > 0).length;
  const open = Control.build(() => button('Open', () => addSavedAnalyses(grok.shell.addTableView(df), analyses),
    {primary: true, tooltip: converted > 0 ?
      `Open the table in the workspace, with ${converted} saved ${converted === 1 ? 'analysis' : 'analyses'}` :
      'Open the table in the workspace'}));
  open.name = 'open';

  return appView({name: fileName, content, ribbon: [[open]]});
}

/** Saved analyses: each one's name, and what it opens as or why it doesn't. Expanded, since unlike the
 * skipped columns this is what the user came to see. */
function analysesSection(analyses: JmpScriptConversion[]): Section {
  const converted = analyses.filter((a) => a.viewers.length > 0).length;
  const section = new Section({
    title: 'Saved analyses',
    hint: `${converted} of ${analyses.length} open as viewers`,
    expanded: true,
    persistKey: 'jmp.preview.analyses',
  });
  section.name = 'saved analyses';
  section.add(divV(analyses.map((a) => {
    const status = a.viewers.length > 0 ? a.viewers.map((v) => v.type).join(', ') : `Not converted: ${a.reason}`;
    const parts = [span(a.name, 'jmp-analysis-name'), span(status, 'jmp-analysis-status')];
    if (a.notes.length > 0)
      parts.push(span(a.notes.join('; '), 'jmp-analysis-notes'));
    const row = divH(parts);
    row.style.gap = '12px';
    row.dataset.u2Item = a.name;
    return row;
  })));
  section.root.style.flex = 'none';
  return section;
}

/** What a file the reader can't open shows instead of a table: why, in the view's own content. */
export function createJmpErrorView(fileName: string, reason: string): DG.ViewBase {
  const content = Control.build(() => {
    const root = divV([
      (() => {
        const title = divH([icon('exclamation-triangle', {size: 'large'}),
          span('This file can\'t be opened as a JMP table', 'jmp-error-title')]);
        Object.assign(title.style, {gap: '8px', fontWeight: 'bold'});
        return title;
      })(),
      span(reason, 'jmp-error-reason'),
      span('JMP data tables saved by current JMP versions are supported. Re-saving the file in JMP may help.',
        'jmp-error-hint'),
    ], 'jmp-error');
    root.dataset.u2Name = 'error';
    Object.assign(root.style, {gap: '8px', padding: '24px', alignItems: 'flex-start'});
    return root;
  });
  return appView({name: fileName, content});
}
