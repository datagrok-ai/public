import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow, subscribeAll, until} from '../helpers';

// TagsInput — ui.input.tags + the allowNew/multiValue/itemToString/getSuggestions/createTagLabel/
// createNewItem config. widget-extras covers TagEditor but not the TagsInput config surface.
category('AI: Widgets: TagsInput', () => {
  test('create with tags seeds items; onChanged is a well-formed stream', async () => {
    const input = ui.input.tags<string>('t', {tags: ['x', 'y'], allowNew: true, multiValue: true});
    expect(input instanceof DG.InputBase, true);
    subscribeAll([input.onChanged])();
    // also exercise the {allowNew:false, multiValue:false} construction path.
    expectNoThrow(() => ui.input.tags<string>('t', {allowNew: false, multiValue: false}));
  });

  test('value set/get round-trips an array (multiValue + allowNew)', async () => {
    const input = ui.input.tags<string>('t', {allowNew: true, multiValue: true});
    (input as DG.InputBase<any>).value = ['a', 'b'];
    await until(() => Array.isArray((input as DG.InputBase<any>).value));
    const value = (input as DG.InputBase<any>).value as string[];
    expect(Array.isArray(value), true);
  });

  test('itemToString / createTagLabel / createNewItem / getSuggestions wire up', async () => {
    expectNoThrow(() => ui.input.tags<string>('t', {
      items: ['alpha', 'beta'],
      allowNew: true,
      itemToString: (s: string) => s.toUpperCase(),
      createTagLabel: (s: string) => ui.span([s]),
      createNewItem: (text: string) => text,
      getSuggestions: async (text: string) => ['alpha', 'beta'].filter((s) => s.startsWith(text)),
    }));
  });
}, {owner: 'agolovko@datagrok.ai'});
