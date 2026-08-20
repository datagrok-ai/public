import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {showProgramDialog, showRejectDialog} from '../app/admin-dialogs';
import {AlignmentHandler, HIDDEN_GRID_COLUMNS} from '../app/catalog-app';
import {T_PROGRAM} from '../domain/constants';
import {cleanupTestPrograms, makeTestProgram} from './fixtures';

category('ArtifactAlignment: ui forms', () => {
  test('program dialog: relevant fields only, derived group names, OK gating', async () => {
    const dialog = showProgramDialog();
    try {
      const captions = dialog.inputs.map((i) => i.caption);
      expect(captions.join(','),
        'Code,Name,Description,External ref,Stage,Status,Indication');
      const ok = dialog.getButton('OK');
      expect(ok.disabled, true, 'OK must be disabled while code and name are empty');
      const code = dialog.input('Code');
      const name = dialog.input('Name');
      code.value = 'ZZUI';
      name.value = 'UI test';
      code.fireChanged();
      name.fireChanged();
      expect(ok.disabled, false, 'OK must enable once code and name are set');
      expect((dialog.root.textContent ?? '').includes('Artifacts: ZZUI Viewers'), true,
        'derived audience group names must be shown');
    } finally {
      dialog.close();
    }
  });

  test('program dialog in edit mode locks the code (business key, group-name seed)', async () => {
    const dialog = showProgramDialog({code: 'ZZED', name: 'Existing'});
    try {
      expect(dialog.input('Code').enabled, false);
      expect(dialog.input('Code').value, 'ZZED');
      expect(dialog.getButton('OK').disabled, false);
    } finally {
      dialog.close();
    }
  });

  test('reject dialog requires a reason', async () => {
    const dialog = showRejectDialog('00000000-0000-0000-0000-000000000000');
    try {
      const ok = dialog.getButton('OK');
      expect(ok.disabled, true, 'OK must be disabled while the reason is empty');
      const reason = dialog.input('Reason');
      reason.value = 'wrong cohort';
      reason.fireChanged();
      expect(ok.disabled, false);
    } finally {
      dialog.close();
    }
  });

  test('program dialog in edit mode links the audience groups', async () => {
    const program = await makeTestProgram();
    const row = await grok.dapi.domains.table(T_PROGRAM).get(program.id);
    const dialog = showProgramDialog(row);
    try {
      await new Promise((resolve) => setTimeout(resolve, 2000));
      const links = Array.from(dialog.root.querySelectorAll('.d4-link-label')).map((e) => e.textContent);
      for (const group of [program.viewers, program.contributors, program.approvers]) {
        expect(links.includes(group.friendlyName), true,
          `group link '${group.friendlyName}' missing, got: ${links.join(', ')}`);
      }
    } finally {
      dialog.close();
    }
  });

  test('catalog grid hides internal bookkeeping columns', async () => {
    const df = DG.DataFrame.fromColumns(
      ['name', ...HIDDEN_GRID_COLUMNS].map((n) => DG.Column.fromStrings(n, ['x'])));
    const grid = DG.Viewer.grid(df);
    new AlignmentHandler().renderGrid(grid);
    for (const name of HIDDEN_GRID_COLUMNS)
      expect(grid.columns.byName(name)?.visible, false, `'${name}' must be hidden`);
    expect(grid.columns.byName('name')?.visible, true, `'name' must stay visible`);
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
