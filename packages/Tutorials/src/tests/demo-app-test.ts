import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect, after, waitUntil} from '@datagrok-libraries/test/src/test';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import type {DemoEntryFunc, DemoView} from '../demo-app/demo-app';

category('Demo App', () => {
  const demoFunc = (name: string, apply: () => Promise<any>): DemoEntryFunc =>
    ({name, description: '', options: {}, package: {name: 'Tutorials'}, apply});
  const openDemoApp = async () => await grok.functions.call('Tutorials:demoApp') as DemoView;
  const addDemoTable = (name: string) => grok.shell.addTableView(DG.DataFrame.create(3, name));

  after(async () => {
    for (const v of Array.from(grok.shell.views).filter((v) => v.temp?.['demoApp'] || v.name === 'user view'))
      v.close();
  });

  test('Demo that opens no view does not rename the user view', async () => {
    const demo = await openDemoApp();
    const userView = grok.shell.addTableView(grok.data.demo.demog(10));
    userView.name = 'user view';
    await demo.startDemoFunc(demoFunc('demoA', async () => addDemoTable('A')), 'Visualization | Demo A');
    await demo.startDemoFunc(demoFunc('demoB', async () => {}), 'Visualization | Demo B');
    expect(userView.name, 'user view');
    expect(userView.path.includes('Demo-B'), false);
  });

  test('Demo selected while another loads shows its own table', async () => {
    const demo = await openDemoApp();
    let releaseX = () => {};
    const xLoading = new Promise<void>((resolve) => releaseX = resolve);
    const x = demo.startDemoFunc(demoFunc('demoX', async () => {
      await xLoading;
      addDemoTable('X');
    }), 'Visualization | Demo X');
    const y = demo.startDemoFunc(demoFunc('demoY', async () => addDemoTable('Y')), 'Visualization | Demo Y');
    releaseX();
    await Promise.all([x, y]);
    expect(grok.shell.tv?.dataFrame.name, 'Y');
    expect(grok.shell.tv?.name, 'Demo Y');
    expect(Array.from(grok.shell.tableViews).some((v) => v.dataFrame.name === 'X'), false);
  });

  test('Demo script next step does not rename the view the user switched to', async () => {
    const userView = grok.shell.addTableView(grok.data.demo.demog(10));
    userView.name = 'user view';
    const otherView = grok.shell.addTableView(grok.data.demo.demog(10));
    otherView.name = 'user view';
    const script = new DemoScript('Demo S', '', false, {autoStartFirstStep: true, path: 'Visualization/Demo S'})
      .step('Open', async () => {addDemoTable('S').temp['demoApp'] = 'demoS';})
      .step('Next', async () => {});
    await script.start();
    const header = document.querySelector('.tutorials-demo-script-container .tutorials-root-header') as HTMLElement;
    try {
      grok.shell.v = userView;
      (header.querySelector('.fa-play')!.closest('button') as HTMLButtonElement).click();
      await waitUntil(() => header.querySelector('.fa-redo') !== null, 'demo script did not finish');
      expect(userView.name, 'user view');
      expect(userView.path.includes('Demo-S'), false);
    }
    finally {
      (header.querySelector('.fa-chevron-left')!.closest('button') as HTMLButtonElement).click();
    }
  });
});
