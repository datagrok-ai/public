import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {subscribeAll, until} from '../helpers';

// TagsInput — ui.input.tags + the allowNew/multiValue/itemToString/getSuggestions/createTagLabel/
// createNewItem config. widget-extras covers TagEditor but not the TagsInput config surface.
category('AI: Widgets: TagsInput', () => {
  test('create with {tags, allowNew, multiValue} builds a working multi-value input', async () => {
    // The {tags} option seeds the available suggestions, not the selected value; the selected value
    // is driven through the setter and is what round-trips.
    const input = ui.input.tags<string>('t', {tags: ['x', 'y'], allowNew: true, multiValue: true});
    expect(input instanceof DG.InputBase, true);
    subscribeAll([input.onChanged])();
    (input as DG.InputBase<any>).value = ['x', 'y'];
    await until(() => Array.isArray((input as DG.InputBase<any>).value) &&
      ((input as DG.InputBase<any>).value as string[]).length === 2);
    const value = (input as DG.InputBase<any>).value as string[];
    expect(value.length, 2);
    expect(value.indexOf('x') >= 0, true);
    expect(value.indexOf('y') >= 0, true);
  });

  test('value set/get round-trips an array (multiValue + allowNew)', async () => {
    const input = ui.input.tags<string>('t', {allowNew: true, multiValue: true});
    (input as DG.InputBase<any>).value = ['a', 'b'];
    await until(() => Array.isArray((input as DG.InputBase<any>).value));
    const value = (input as DG.InputBase<any>).value as string[];
    expect(Array.isArray(value), true);
    expect(value.length, 2);
    expect(value[0], 'a');
    expect(value[1], 'b');
  });

  test('itemToString / createTagLabel / createNewItem / getSuggestions seed items and round-trip value', async () => {
    const input = ui.input.tags<string>('t', {
      items: ['alpha', 'beta'],
      allowNew: true,
      multiValue: true,
      itemToString: (s: string) => (s ?? '').toUpperCase(),
      createTagLabel: (s: string) => ui.div([s]),
      createNewItem: (text: string) => text,
      getSuggestions: async (text: string) => ['alpha', 'beta'].filter((s) => s.startsWith(text)),
    });
    (input as DG.InputBase<any>).value = ['alpha', 'beta'];
    await until(() => Array.isArray((input as DG.InputBase<any>).value));
    const value = (input as DG.InputBase<any>).value as string[];
    expect(value.length, 2);
    expect(value.indexOf('alpha') >= 0, true);
    expect(value.indexOf('beta') >= 0, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
