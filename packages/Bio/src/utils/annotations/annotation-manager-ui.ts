/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  SeqAnnotation, AnnotationCategory, LiabilitySeverity,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {getColumnAnnotations, setColumnAnnotations, clearAnnotations} from './annotation-manager';

const categoryLabels: Record<string, string> = {
  [AnnotationCategory.Structure]: 'Structure (FR/CDR)',
  [AnnotationCategory.Liability]: 'Liability',
  [AnnotationCategory.PTM]: 'Post-translational Modification',
  [AnnotationCategory.Custom]: 'Custom',
};

const severityLabels: Record<string, string> = {
  [LiabilitySeverity.High]: 'High',
  [LiabilitySeverity.Medium]: 'Medium',
  [LiabilitySeverity.Low]: 'Low',
  [LiabilitySeverity.Info]: 'Info',
};

export function showAnnotationManagerDialog(): void {
  const df = grok.shell.tv?.dataFrame;
  if (!df) {
    grok.shell.warning('No table open');
    return;
  }

  const seqCols = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (seqCols.length === 0) {
    grok.shell.warning('No macromolecule columns found');
    return;
  }

  let selectedCol = seqCols[0];
  const colInput = ui.input.column('Sequence Column', {
    table: df, value: selectedCol,
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
    onValueChanged: (col) => { selectedCol = col!; refreshList(); },
  });

  const listDiv = ui.divV([], {style: {maxHeight: '380px', overflowY: 'auto', paddingRight: '8px'}});

  function refreshList(): void {
    listDiv.innerHTML = '';
    const annotations = getColumnAnnotations(selectedCol);
    if (annotations.length === 0) {
      listDiv.append(ui.divText('No annotations on this column.', {style: {color: '#888', padding: '8px'}}));
      return;
    }

    for (const annot of annotations) {
      const catLabel = categoryLabels[annot.category] ?? annot.category;
      const sevLabel = annot.severity ? ` [${severityLabels[annot.severity] ?? annot.severity}]` : '';
      const rangeLabel = annot.start && annot.end ? ` (${annot.start}-${annot.end})` : '';
      const schemeLabel = annot.sourceScheme ? ` ${annot.sourceScheme}` : '';

      const removeBtn = ui.iconFA('trash', () => {
        const updated = getColumnAnnotations(selectedCol).filter((a) => a.id !== annot.id);
        setColumnAnnotations(selectedCol, updated);
        df.fireValuesChanged();
        refreshList();
      });
      removeBtn.style.cursor = 'pointer';
      removeBtn.style.color = '#999';
      removeBtn.style.marginLeft = '8px';
      const originalColor = annot.color ?? '#ccc';
      let currentColor = originalColor;
      const colorSwatch = ui.div([], {style: {
        width: '12px', height: '12px', borderRadius: '2px',
        backgroundColor: currentColor, display: 'inline-block', marginRight: '6px',
        flexShrink: '0', cursor: 'pointer',
      }});

      ui.colorPicker(DG.Color.fromHtml(annot.color ?? '#ccc'), (newColor) => {
        currentColor = DG.Color.toHtml(newColor);
      }, colorSwatch, () => {
        const updated = getColumnAnnotations(selectedCol).map((a) => a.id === annot.id ? {...a, color: currentColor} : a);
        setColumnAnnotations(selectedCol, updated);
        df.fireValuesChanged();
        refreshList();
      }, () => {
        currentColor = originalColor;
        colorSwatch.style.backgroundColor = currentColor;
      });

      const row = ui.divH([
        colorSwatch,
        ui.divText(`${annot.name}${rangeLabel}${schemeLabel}${sevLabel}`, {style: {flex: '1', fontSize: '12px', padding: '4px'}}),
        ui.divText(catLabel, {style: {color: '#888', fontSize: '11px', marginRight: '8px'}}),
        removeBtn,
      ], {style: {alignItems: 'center', padding: '4px 0', borderBottom: '1px solid #eee'}});

      listDiv.append(row);
    }
  }

  refreshList();

  const clearBtn = ui.button('Clear All', () => {
    clearAnnotations(df, selectedCol);
    df.fireValuesChanged();
    refreshList();
    grok.shell.info('All annotations cleared');
  });

  const dialog = ui.dialog({title: 'Manage Annotations'})
    .add(ui.inputs([colInput]))
    .add(ui.h3('Annotations'))
    .add(listDiv)
    .add(ui.divH([clearBtn], {style: {marginTop: '8px'}}))
    .onOK(() => {});

  dialog.show();
}
