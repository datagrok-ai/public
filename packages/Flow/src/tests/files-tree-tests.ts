/** The Files tree (KNIME-style browser) stamps a stable `data-testid` on every
 *  connection, folder, and file row — the hook the guides use to steer a user to
 *  drag a real file (e.g. Demo / demog.csv) onto the canvas. Needs a live stand
 *  (lists connections + files). */
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {getFilesBrowser} from '../utils/files-browser-tree';
import {FuncFlowView} from '../funcflow-view';
import {ensureFunctionsRegistered} from '../rete/node-factory';

const sleep = (ms: number): Promise<void> => new Promise((r) => setTimeout(r, ms));

async function waitFor<T>(fn: () => T | null | undefined, timeoutMs = 8000): Promise<T | null> {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    const v = fn();
    if (v) return v;
    await sleep(120);
  }
  return null;
}

category('Flow: files tree', () => {
  test('a dropped file titles its Open File node with the file name', async () => {
    // Several Open File nodes on one canvas must be tellable apart at a glance
    // — the title carries the file name (cosmetic; never affects the run).
    ensureFunctionsRegistered();
    const view = new FuncFlowView();
    const probe = view as unknown as {
      editorReady: Promise<void>;
      flow: {getNodes(): Array<{label: string; inputValues: Record<string, unknown>}>; destroy?(): void};
      addOpenFileNode(path: string): Promise<void>;
    };
    try {
      await probe.editorReady;
      await probe.addOpenFileNode('System:AppData/Chem/mol1K.csv');
      const node = probe.flow.getNodes()[0];
      expect(node != null, true, 'node added');
      expect(node.label.endsWith(': mol1K.csv'), true, `title carries the file name (got "${node.label}")`);
      expect(node.inputValues['fullPath'], 'System:AppData/Chem/mol1K.csv', 'path still configured');
    } finally {
      // The path write arms the live-by-default autorun (1 s debounce) —
      // detach resets the scheduler so nothing fires into the torn-down view.
      view.detach();
      await sleep(120);
      probe.flow?.destroy?.();
      view.root.remove();
    }
  }, {timeout: 30000});

  test('connections, folders and files carry name-based test-ids', async () => {
    const tree = getFilesBrowser(() => {}, () => {}, undefined);
    document.body.appendChild(tree.root);
    try {
      // Connections load asynchronously; the Demo connection is always present.
      const demoGroup = await waitFor(() =>
        tree.children.find((c): c is DG.TreeViewGroup =>
          c instanceof DG.TreeViewGroup && (c.root as HTMLElement).dataset.conn === 'Demo') ?? null);
      expect(!!demoGroup, true, 'Demo connection row present');
      const demoRoot = demoGroup!.root as HTMLElement;
      expect(demoRoot.dataset.testid, 'ff-files-conn-demo', 'connection test-id is name-based');

      // Expanding the connection lazily loads + stamps its files.
      demoGroup!.expanded = true;
      const fileRow = await waitFor(() =>
        tree.root.querySelector('[data-testid="ff-files-file-demog-csv"]') as HTMLElement | null);
      expect(!!fileRow, true, 'demog.csv file row present with name-based test-id');
      expect((fileRow!.dataset.file ?? '').toLowerCase(), 'demog.csv', 'raw file name attribute');
      expect((fileRow!.dataset.filePath ?? '').toLowerCase().includes('demog'), true, 'fullPath attribute set');

      // A folder row is stamped too (e.g. the Demo/chem folder).
      const folderRow = tree.root.querySelector('[data-testid="ff-files-folder-chem"]') as HTMLElement | null;
      expect(!!folderRow, true, 'folder row stamped with a name-based test-id');
    } finally {
      tree.root.remove();
    }
  });
});
