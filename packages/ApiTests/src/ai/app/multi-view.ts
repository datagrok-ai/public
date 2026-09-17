import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, uniqueName, wait} from '../helpers';

category('AI: App: MultiView', () => {
  // MultiView is a ViewBase, not a View, so withAttachedView() does not accept it.
  async function withMultiView(factories: {[name: string]: () => DG.ViewBase},
    body: (mv: DG.MultiView) => Promise<void>): Promise<void> {
    const mv = new DG.MultiView({viewFactories: factories});
    grok.shell.addView(mv);
    try {
      await body(mv);
    } finally {
      mv.close();
    }
  }

  const activate = async (mv: DG.MultiView, name: string, ms: number = 300): Promise<void> => {
    mv.tabs.currentPane = mv.tabs.getPane(name);
    await wait(ms);
  };

  // Handing a ribbon item to another view moves it out of its `d4-ribbon-item` wrapper; reading
  // the ribbon again used to return the empty wrappers, so the buttons vanished on the way back.
  test('ribbon buttons survive a tab round-trip', async () => {
    const button = ui.button('hey', () => {});
    const plain = DG.View.fromRoot(ui.divText('custom content'));
    plain.setRibbonPanels([[button]]);

    await withMultiView({
      'plain': () => plain,
      'other': () => DG.View.fromRoot(ui.divText('other')),
    }, async (mv) => {
      const shown = (): boolean => mv.root.parentElement!.querySelector('.d4-ribbon')!.contains(button);

      await activate(mv, 'plain');
      expect(shown(), true);

      await activate(mv, 'other');
      await activate(mv, 'plain');
      expect(shown(), true);
    });
  });

  // shell.t resolves through ViewBase.innerView, so commands gated on a current table view
  // (Select | All, Add New Column, ...) act on the view inside the tab rather than being disabled.
  test('shell.t resolves into the active tab', async () => {
    const df = demog();
    df.name = uniqueName('mv-table');

    await withMultiView({
      'table': () => DG.TableView.create(df, false),
      'plain': () => DG.View.fromRoot(ui.divText('no table here')),
    }, async (mv) => {
      await activate(mv, 'table');
      expect(grok.shell.t?.name, df.name);

      await activate(mv, 'plain');
      expect(grok.shell.t == null, true);
    });
  });

  // A Dart view root carries no ui layout class until the shell opens it, so hosted in a tab it
  // collapsed to the height of its header and the grid never rendered.
  test('a hosted table view fills its tab', async () => {
    await withMultiView({
      'table': () => DG.TableView.create(demog(), false),
    }, async (mv) => {
      await activate(mv, 'table', 600);
      const pane = mv.currentView.root.parentElement!;
      expect(mv.currentView.root.clientHeight > pane.clientHeight / 2, true);
    });
  });
});
