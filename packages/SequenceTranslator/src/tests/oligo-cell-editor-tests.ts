import {after, category, test, expect, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {tagAsOligoNucleotide} from '../oligo-renderer/converters';

const SAMPLE_HELM =
  'RNA1{r(A)p.r(C)p.r(G)p.r(U)p}|RNA2{r(U)p.r(C)p.r(G)p.r(A)p}$$$$';

function dialogCount(): number {
  return $('.d4-dialog').length;
}

function closeAllDialogs(): void {
  $('.d4-dialog .ui-btn-cancel, .d4-dialog .d4-dialog-header .grok-icon.fa-times').trigger('click');
}

category('OligoCellEditor', () => {
  after(async () => {
    closeAllDialogs();
    grok.shell.closeAll();
  });

  test('cellEditor opens HELM editor for OligoNucleotide cell and saves on OK', async () => {
    const col = DG.Column.fromStrings('seq', [SAMPLE_HELM]);
    tagAsOligoNucleotide(col);
    const df = DG.DataFrame.fromColumns([col]);
    df.name = 'oligo-edit-test';
    const tv = grok.shell.addTableView(df);

    await awaitCheck(() => $(tv.root).find('.d4-grid canvas').length > 0,
      'Grid canvas did not appear', 5000);

    // Find the cellEditor that the platform would dispatch on double-click for
    // a column tagged quality=OligoNucleotide. Pre-fix: zero matches (the
    // registration didn't exist). Post-fix: exactly one — editOligoNucleotideCell.
    const matches = DG.Func.find({tags: ['cellEditor'], package: 'SequenceTranslator'})
      .filter((f) => f.description === 'OligoNucleotide');
    expect(matches.length, 1);
    expect(matches[0].name, 'editOligoNucleotideCell');

    const gridCell = tv.grid.cell('seq', 0);
    expect(gridCell != null, true);
    expect(gridCell.cell.value, SAMPLE_HELM);

    const dialogsBefore = dialogCount();

    // Invoke the cellEditor as the platform would on double-click.
    // Pre-fix (delegating to Helm:editMoleculeCell): throws synchronously
    //   "The column of notation 'helm' must be 'Macromolecule'" — dialog never opens.
    // Post-fix (using helmHelper.createWebEditorApp directly): dialog opens.
    await matches[0].apply({cell: gridCell});

    await awaitCheck(() => dialogCount() > dialogsBefore,
      'HELM editor dialog did not open within 15s', 15000);

    // Wait for HWE async init (Dojo + JSDraw2 + monomer lib) to mount the editor.
    // JSDraw2 renders to SVG, so wait for the OK button to be wired up — that's
    // a reliable signal that the dialog footer is fully constructed.
    await awaitCheck(() => $('.d4-dialog .ui-btn-ok, .d4-dialog button.ui-btn.ui-btn-ok').length > 0,
      'OK button did not appear in HELM editor dialog within 15s', 15000);

    // Allow the editor a moment to load the HELM string into the canvas before we read it back.
    await delay(1000);

    const okBtn = $('.d4-dialog .ui-btn-ok, .d4-dialog button.ui-btn.ui-btn-ok').first();
    expect(okBtn.length > 0, true, 'OK button not found in dialog');

    okBtn.trigger('click');

    await awaitCheck(() => dialogCount() <= dialogsBefore,
      'Dialog did not close after OK', 5000);

    // After OK: cell.setValue(helmValue) must have run with the editor's HELM.
    // The editor may canonicalize formatting, so we don't require byte-equality —
    // we require it remains a valid two-strand HELM string for our input.
    const after = gridCell.cell.value as string;
    expect(typeof after === 'string' && after.includes('RNA1{') && after.includes('RNA2{') && after.includes('$$$$'), true,
      `Expected cell value to remain valid HELM after OK; got: ${after}`);
  });
});
