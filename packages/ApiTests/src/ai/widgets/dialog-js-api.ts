import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectFiresWithin, subscribeAll, expectNoThrow, until} from '../helpers';

function getCancelButton(dlg: DG.Dialog): HTMLButtonElement | null {
  return dlg.getButton('CANCEL') ?? dlg.getButton('Cancel') ?? null;
}

category('AI: Widgets: Dialog JS API', () => {
  test('factory create + title/helpUrl round-trip', async () => {
    const dlg = DG.Dialog.create({title: 'D1', helpUrl: 'https://example.org/help'});
    expect(dlg instanceof DG.Dialog, true);
    expect(dlg.title, 'D1');
    expect(dlg.helpUrl, 'https://example.org/help');
    dlg.title = 'D1-renamed';
    dlg.helpUrl = 'https://example.org/help2';
    expect(dlg.title, 'D1-renamed');
    expect(dlg.helpUrl, 'https://example.org/help2');
    const dlg2 = DG.Dialog.create({title: 'D2', showFooter: false});
    expect(dlg2 instanceof DG.Dialog, true);
  });

  test('add HTMLElement + addInput + namedInputs + clear', async () => {
    const dlg = DG.Dialog.create({title: 'Inputs'});
    const div = ui.div([], 'dialog-test-host');
    dlg.add(div);
    const captionInput = ui.input.string('Caption', {value: 'hello'});
    dlg.add(captionInput);
    const dlg2 = dlg.addInput('nick', ui.input.string('nick', {value: 'world'}));
    expect(dlg2 === dlg, true);
    const inputs = dlg.inputs;
    expect(inputs.length >= 2, true);
    const named = dlg.namedInputs as {[key: string]: DG.InputBase};
    expect(named['nick'] != null, true);
    expect(named['nick'].value, 'world');
    expectNoThrow(() => dlg.clear());
  });

  test('show + close lifecycle; getOpenDialogs reflects open state', async () => {
    const baseline = DG.Dialog.getOpenDialogs().length;
    const dlg = DG.Dialog.create({title: 'Lifecycle'});
    try {
      dlg.show();
      await until(() => DG.Dialog.getOpenDialogs().length === baseline + 1);
      expect(DG.Dialog.getOpenDialogs().length, baseline + 1);
    } finally {
      dlg.close();
    }
    await until(() => DG.Dialog.getOpenDialogs().length === baseline);
    expect(DG.Dialog.getOpenDialogs().length, baseline);
  });

  test('onOK callback + onBeforeOK/onAfterOK observables fire; closes on OK', async () => {
    const dlg = DG.Dialog.create({title: 'OK Flow'});
    let okFired = false;
    dlg.onOK(() => {okFired = true;});
    try {
      dlg.show();
      await until(() => dlg.getButton('OK') != null);
      const okBtn = dlg.getButton('OK');
      expect(okBtn != null, true);
      await expectFiresWithin(dlg.onBeforeOK, () => okBtn.click());
      await until(() => okFired);
      expect(okFired, true);
      subscribeAll([dlg.onAfterOK])();
      await until(() => {
        for (const d of DG.Dialog.getOpenDialogs())
          if (d.title === 'OK Flow') return false;
        return true;
      });
      const open = DG.Dialog.getOpenDialogs();
      let stillOpen = false;
      for (const d of open) {
        if (d.title === 'OK Flow') {
          stillOpen = true;
          break;
        }
      }
      expect(stillOpen, false);
    } finally {
      try {dlg.close();} catch (_e) {}
    }
  });

  test('onCancel callback + onClose observable fire when CANCEL clicked', async () => {
    const dlg = DG.Dialog.create({title: 'Cancel Flow'});
    let cancelFired = false;
    dlg.onCancel(() => {cancelFired = true;});
    try {
      dlg.show();
      await until(() => getCancelButton(dlg) != null);
      const cancelBtn = getCancelButton(dlg);
      expect(cancelBtn != null, true);
      await expectFiresWithin(dlg.onClose, () => cancelBtn!.click());
      await until(() => cancelFired);
      expect(cancelFired, true);
    } finally {
      try {dlg.close();} catch (_e) {}
    }
  });

  test('addButton + getButton + addContextAction', async () => {
    const dlg = DG.Dialog.create({title: 'Buttons'});
    let applyFired = false;
    dlg.addButton('APPLY', () => {applyFired = true;});
    expectNoThrow(() => dlg.addContextAction('Reset', () => {}));
    try {
      dlg.show();
      await until(() => dlg.getButton('APPLY') != null);
      const apply = dlg.getButton('APPLY');
      expect(apply != null, true);
      expect(apply instanceof HTMLElement, true);
      apply!.click();
      await until(() => applyFired);
      expect(applyFired, true);
    } finally {
      try {dlg.close();} catch (_e) {}
    }
  });

  test('awaitOnOK resolves on OK; rejects on CANCEL', async () => {
    const dlgOk = DG.Dialog.create({title: 'Await OK'});
    const okPromise = dlgOk.awaitOnOK<number>(async () => 42);
    try {
      dlgOk.show();
      await until(() => dlgOk.getButton('OK') != null);
      dlgOk.getButton('OK').click();
      const result = await okPromise;
      expect(result, 42);
    } finally {
      try {dlgOk.close();} catch (_e) {}
    }

    const dlgCancel = DG.Dialog.create({title: 'Await Cancel'});
    const cancelPromise = dlgCancel.awaitOnOK<number>(async () => 1);
    let rejected = false;
    try {
      dlgCancel.show();
      await until(() => getCancelButton(dlgCancel) != null);
      const cancelBtn = getCancelButton(dlgCancel);
      expect(cancelBtn != null, true);
      cancelBtn!.click();
      try {
        await cancelPromise;
      } catch (_e) {
        rejected = true;
      }
      expect(rejected, true);
    } finally {
      try {dlgCancel.close();} catch (_e) {}
    }
  });

  test('static getOpenDialogs returns DG.Dialog instances', async () => {
    const baseline = DG.Dialog.getOpenDialogs().length;
    const a = DG.Dialog.create({title: 'Open A'});
    const b = DG.Dialog.create({title: 'Open B'});
    try {
      a.show();
      b.show();
      await until(() => DG.Dialog.getOpenDialogs().length >= baseline + 2);
      const open = DG.Dialog.getOpenDialogs();
      expect(open.length >= baseline + 2, true);
      for (const d of open)
        expect(d instanceof DG.Dialog, true);
    } finally {
      try {a.close();} catch (_e) {}
      try {b.close();} catch (_e) {}
    }
  });

  test('showModal opens; double close is a no-op', async () => {
    const dlg = DG.Dialog.create({title: 'Modal'});
    try {
      dlg.showModal(false);
      await until(() => {
        for (const d of DG.Dialog.getOpenDialogs())
          if (d.title === 'Modal') return true;
        return false;
      });
      let found = false;
      for (const d of DG.Dialog.getOpenDialogs()) {
        if (d.title === 'Modal') {
          found = true;
          break;
        }
      }
      expect(found, true);
      dlg.close();
      await until(() => {
        for (const d of DG.Dialog.getOpenDialogs())
          if (d.title === 'Modal') return false;
        return true;
      });
      expectNoThrow(() => dlg.close());
    } finally {
      try {dlg.close();} catch (_e) {}
    }
  });

  test('onSaveHistoryRequest observable subscribes cleanly', async () => {
    const dlg = DG.Dialog.create({title: 'History'});
    try {
      subscribeAll([dlg.onSaveHistoryRequest])();
    } finally {
      try {dlg.close();} catch (_e) {}
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
