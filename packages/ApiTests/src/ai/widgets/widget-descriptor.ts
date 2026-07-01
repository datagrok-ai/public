import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// DG.WidgetDescriptor — the viewer/widget registry: getDescriptors/getByName,
// name/synonyms/description/properties/events getters, createIcon().
category('AI: Widgets: WidgetDescriptor', () => {
  test('getDescriptors returns a non-empty registry of WidgetDescriptor', async () => {
    const all = DG.WidgetDescriptor.getDescriptors();
    expect(Array.isArray(all), true);
    expect(all.length > 0, true);
    for (const d of all) {
      expect(d instanceof DG.WidgetDescriptor, true);
      expect(typeof d.name, 'string');
      expect(d.name.length > 0, true);
    }
  });

  test('getByName round-trips a known descriptor name', async () => {
    const first = DG.WidgetDescriptor.getDescriptors()[0];
    const found = DG.WidgetDescriptor.getByName(first.name);
    expect(found != null, true);
    expect(found!.name, first.name);
  });

  test('getByName returns null/undefined for an unknown name', async () => {
    expect(DG.WidgetDescriptor.getByName('No Such Widget 12345') == null, true);
  });

  test('descriptor exposes synonyms / description / properties / events', async () => {
    const d = DG.WidgetDescriptor.getDescriptors()[0];
    expect(Array.isArray(d.synonyms), true);
    expect(typeof d.description, 'string');
    expect(Array.isArray(d.properties), true);
    expect(Array.isArray(d.events), true);
  });

  test('createIcon returns a DOM Element', async () => {
    const d = DG.WidgetDescriptor.getDescriptors()[0];
    const icon = d.createIcon();
    expect(icon instanceof Element, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
