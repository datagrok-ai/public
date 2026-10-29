import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {filter, take} from 'rxjs/operators';

import {badge, Control, div, divH, divV, icon, link, Section, show, signal, span} from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';

import {JmpTable} from './jmp-types';
import {addSavedAnalyses, dataColumns, hideJmpHiddenColumns, jmpTableToDataFrame, savedAnalyses, skippedColumns}
  from './jmp-to-dataframe';
import {JmpScriptConversion} from './jsl-to-viewers';

const list = (names: string[]): string => `${names.slice(0, 5).join(', ')}${names.length > 5 ? ', …' : ''}`;

/** Says which columns and saved analyses didn't come across, so nothing disappears silently. */
export function warnSkipped(table: JmpTable): void {
  const names = skippedColumns(table).map((s) => s.name);
  if (names.length > 0)
    grok.shell.warning(`JMP: ${names.length} column${names.length === 1 ? '' : 's'} not imported: ${list(names)}`);
  const analyses = savedAnalyses(table).filter((a) => a.viewers.length === 0).map((a) => a.name);
  if (analyses.length > 0) {
    grok.shell.warning(`JMP: ${analyses.length} saved ${analyses.length === 1 ? 'analysis' : 'analyses'} ` +
      `not converted: ${list(analyses)}`);
  }
}

/** A JMP table as a table view: its data, the saved analyses that convert open as viewers, and a ribbon
 * button for the file summary. Both Browse (click to preview, double-click to keep) and a file opened from
 * disk show this, so what the user previews is what they open (JMP_PLUGIN.md §6.1). */
export function createJmpTableView(table: JmpTable, fileName: string): DG.TableView {
  const tv = DG.TableView.create(jmpTableToDataFrame(table, fileName.replace(/\.jmp$/i, '')), false);
  // A table view has no grid, and takes no viewers, until the shell shows it, as a preview or a view. The
  // context panel's small preview of the file builds a view it never shows, so it doesn't warn twice.
  const shown = grok.events.onViewAdded.pipe(filter((v) => v.dart === tv.dart), take(1)).subscribe(() => {
    warnSkipped(table);
    hideJmpHiddenColumns(tv);
    addSavedAnalyses(tv, savedAnalyses(table));
    const about = ui.iconFA('info-circle', () => showJmpSummary(table, fileName), 'About this JMP file');
    about.dataset.jmpName = 'summary button';
    about.setAttribute('aria-label', 'About this JMP file');
    tv.setRibbonPanels([...tv.getRibbonPanels(), [about]]);
  });
  tv.subs.push(shown);
  return tv;
}

/** What the file is, the columns that didn't come across, and each saved analysis with what it opened as
 * or why it didn't. */
export function showJmpSummary(table: JmpTable, fileName: string): DG.Dialog {
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
      parts.push(section);
    }
    if (analyses.length > 0)
      parts.push(analysesSection(analyses));
    const root = divV(parts, 'jmp-preview');
    Object.assign(root.style, {maxWidth: '800px', maxHeight: '70vh', overflow: 'auto'});
    return root;
  });
  const dialog = ui.dialog(fileName).add(content.root);
  dialog.onClose.subscribe(() => content.dispose());
  return dialog.show();
}

/** Saved analyses: each one's name, what it opens as or why it doesn't, and its script on demand.
 * Expanded, since unlike the skipped columns this is what the user came to see. */
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
    const scriptShown = signal(false);
    const toggle = link('Script', () => scriptShown.value = !scriptShown.value);
    toggle.dataset.u2Name = 'script toggle';
    const parts = [span(a.name, 'jmp-analysis-name'), span(status, 'jmp-analysis-status')];
    if (a.notes.length > 0)
      parts.push(span(a.notes.join('; '), 'jmp-analysis-notes'));
    parts.push(toggle);
    const line = divH(parts);
    line.style.gap = '12px';
    const scriptHost = div();
    show(scriptHost, scriptShown, () => {
      const text = div([a.text], 'jmp-analysis-script');
      text.dataset.u2Name = 'script';
      Object.assign(text.style, {whiteSpace: 'pre-wrap', fontFamily: 'monospace', fontSize: '11px',
        maxHeight: '160px', overflow: 'auto', padding: '4px 8px', userSelect: 'text'});
      return text;
    });
    const row = divV([line, scriptHost]);
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
