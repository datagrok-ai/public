import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, subscribeAll} from '../helpers';

category('AI: Widgets: Extras', () => {
  test('RangeSlider: create + setValues + scrollTo + scrollBy', async () => {
    const rs = DG.RangeSlider.create();
    expect(rs instanceof DG.RangeSlider, true);
    rs.setValues(0, 100, 10, 90);
    expectNoThrow(() => {
      rs.scrollTo(50);
      rs.scrollBy(5);
    });
  });

  test('ColumnComboBox: create + vertical + onChanged + onEvent', async () => {
    const ccb = DG.ColumnComboBox.create(demog(), (c: DG.Column) => c.type === 'double');
    expect(ccb instanceof DG.ColumnComboBox, true);
    expect(typeof ccb.vertical, 'boolean');
    const evt = ccb.onEvent('d4-column-box-column-changed');
    expect(typeof evt.subscribe, 'function');
    subscribeAll([ccb.onChanged])();
  });

  test('Legend: create + selected categories + tooltip mode + onViewerLegendChanged', async () => {
    const legend = DG.Legend.create(demog().col('sex')!);
    expect(legend instanceof DG.Legend, true);
    expectNoThrow(() => {
      const selected = legend.selectedCategories;
      expect(selected == null || Array.isArray(selected), true);
      const extra = legend.selectedExtraCategories;
      expect(extra == null || Array.isArray(extra), true);
      expect(legend.isTooltipMode, false);
      legend.onViewerLegendChanged = () => {};
    });
  });

  test('InputForm: forInputs + isValid + onValidationCompleted', async () => {
    const form = DG.InputForm.forInputs([ui.input.string('a'), ui.input.int('b')]);
    expect(form instanceof DG.InputForm, true);
    expect(typeof form.isValid, 'boolean');
    subscribeAll([form.onValidationCompleted])();
  });

  test('DropDown: menu + isExpanded + expand + collapse + removeMenuSubscription', async () => {
    const dd = DG.DropDown.menu('Actions', {A: () => {}, B: () => {}});
    expect(dd instanceof DG.DropDown, true);
    expect(dd.isExpanded, false);
    dd.expand();
    expect(dd.isExpanded, true);
    dd.collapse();
    expect(dd.isExpanded, false);
    expectNoThrow(() => dd.removeMenuSubscription());
  });

  test('TagEditor: create + addTag/removeTag/clearTags + drag-drop setters + onChanged', async () => {
    const te = DG.TagEditor.create();
    expect(te instanceof DG.TagEditor, true);
    const sub = te.onChanged(() => {});
    try {
      te.addTag('x');
      te.addTag('y');
      te.removeTag('x');
      te.clearTags();
      expectNoThrow(() => {
        te.acceptsDragDrop = (..._args: any[]) => true;
        te.doDrop = (..._args: any[]) => {};
      });
    } finally {
      sub.cancel();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
