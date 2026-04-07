import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {
  BUILTIN_LIABILITY_RULES, LiabilityRule, scanLiabilities,
  applyLiabilityScanResults, createLiabilitySummaryColumn,
} from './liability-scanner';
import {LiabilitySeverity} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';

const severityLabels: Record<string, string> = {
  [LiabilitySeverity.High]: 'High',
  [LiabilitySeverity.Medium]: 'Medium',
  [LiabilitySeverity.Low]: 'Low',
  [LiabilitySeverity.Info]: 'Info',
};

export function showLiabilityScannerDialog(): void {
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

  const rules = BUILTIN_LIABILITY_RULES.map((r) => ({...r, pattern: new RegExp(r.pattern.source, 'g')}));

  const tableInput = ui.input.table('Table', {value: df});
  const seqInput = ui.input.column('Sequence', {
    table: df, value: seqCols[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
  });

  // Rule checkboxes
  const ruleChecks: {rule: LiabilityRule; input: DG.InputBase<boolean>}[] = [];
  const rulesDiv = ui.divV([]);
  for (const rule of rules) {
    const check = ui.input.bool(rule.name, {
      value: rule.enabled,
      tooltipText: `Severity: ${severityLabels[rule.severity] ?? rule.severity}`,
    });
    ruleChecks.push({rule, input: check});
    rulesDiv.append(check.root);
  }

  const highlightInput = ui.input.bool('Highlight in cell renderer', {value: true});
  const annotColInput = ui.input.bool('Create annotation column', {value: true});
  const summaryInput = ui.input.bool('Create summary count column', {value: false});

  const dialog = ui.dialog({title: 'Scan Sequence Liabilities'})
    .add(ui.inputs([tableInput, seqInput]))
    .add(ui.h3('Rules'))
    .add(rulesDiv)
    .add(ui.h3('Output'))
    .add(ui.inputs([highlightInput, annotColInput, summaryInput]))
    .onOK(() => {
      try {
        const seqCol = seqInput.value!;
        const sh = _package.seqHelper.getSeqHandler(seqCol);

        // Apply checkbox state
        for (const {rule, input} of ruleChecks)
          rule.enabled = input.value;

        const result = scanLiabilities(seqCol, sh, rules);

        if (annotColInput.value || highlightInput.value)
          applyLiabilityScanResults(df, seqCol, result);

        if (summaryInput.value)
          createLiabilitySummaryColumn(df, seqCol, result);

        grok.shell.info(`Liability scan: ${result.totalHits} hits found across ${result.annotations.length} rules`);
        df.fireValuesChanged();
      } catch (err: any) {
        grok.shell.error(`Liability scan failed: ${err.message ?? err}`);
        console.error(err);
      }
    });

  dialog.show();
}
