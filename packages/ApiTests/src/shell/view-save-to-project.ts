import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

const nodeSkip = typeof process !== 'undefined' ? 'NodeJS environment' : undefined;

const inScratchpad = (v: DG.ViewBase): boolean =>
  grok.shell.project.children.some((c) => c instanceof DG.ViewInfo && c.name === v.name);

category('View: Save to project', () => {
  test('a markup view waits to be included, and the flag includes it or takes it out', async () => {
    await grok.functions.call('CmdNewMarkupView');
    const markup = grok.shell.v;
    try {
      expect(markup.saveToProjectByDefault, false);
      expect(markup.saveToProject, false);
      expect(inScratchpad(markup), true, 'a new view is offered in the Save project dialog');
      markup.saveToProject = true;
      expect(markup.saveToProject, true);
      expect(inScratchpad(markup), true);
      markup.saveToProject = false;
      expect(markup.saveToProject, false);
      expect(inScratchpad(markup), false);
    } finally {
      markup.close();
    }
  }, {skipReason: nodeSkip});

  test('a view that cannot be reopened is never saved', async () => {
    const view = grok.shell.addView(DG.View.create());
    try {
      view.saveToProject = true;
      expect(view.saveToProject, false);
      expect(inScratchpad(view), false);
    } finally {
      view.close();
    }
  }, {skipReason: nodeSkip});
}, {owner: 'dkovalyov@datagrok.ai'});
