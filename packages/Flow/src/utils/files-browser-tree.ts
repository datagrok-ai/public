/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {setTid} from './test-ids';

type FileArgFunc = (f: DG.FileInfo) => void;

// Stable, addressable test-ids on every tree row so tutorials and UI tests can
// target a specific connection / folder / file by its name (e.g. the Demo
// connection's expander, or the demog.csv file).
const stampConnection = (root: HTMLElement, name: string): void => {
  setTid(root, 'files-conn', name);
  root.dataset.conn = name;
};
const stampFolder = (root: HTMLElement, name: string): void => {
  setTid(root, 'files-folder', name);
  root.dataset.folder = name;
};
const stampFile = (root: HTMLElement, f: DG.FileInfo): void => {
  setTid(root, 'files-file', f.fileName);
  root.dataset.file = f.fileName;
  root.dataset.filePath = f.fullPath;
};

/** Hierarchical snapshot of which groups are expanded, keyed by node key (see `nodeKey`). */
interface ExpandedState {
  [key: string]: ExpandedState;
}

// onSelected will be called when a file node is selected.
// onDblClick will be called when a file node is double-clicked (or activated with Enter).
// An optional localStorage key will be used to save and restore the expanded state of the tree.
export function getFilesBrowser(onSelected: FileArgFunc, onDblClick: FileArgFunc, localStorageKey?: string): DG.TreeViewGroup {
  const tree = ui.tree();

  const addLoader = (g: DG.TreeViewGroup) => {
    const loaderElem: HTMLElement = ui.loader();
    loaderElem.style.width = '50px';
    const item = g.item(loaderElem);
    return () => item.remove();
  };

  // A stable identity for a group, used both as the localStorage key and to match
  // saved state against live nodes. Values are unique within their parent.
  const nodeKey = (node: DG.TreeViewNode): string | null => {
    const v = node.value;
    if (v instanceof DG.FileInfo)
      return v.fileName;
    if (v instanceof DG.DataConnection)
      return v.nqName;
    return null;
  };

  // Lists the content of a node's value: the connection root for a connection node,
  // or the directory itself for a folder node.
  const listContent = async (value: unknown): Promise<DG.FileInfo[]> => {
    if (value instanceof DG.DataConnection) {
      const path = value.nqName.endsWith('/') ? value.nqName : value.nqName + '/';
      return grok.dapi.files.list(path);
    }
    if (value instanceof DG.FileInfo && value.isDirectory)
      return grok.dapi.files.list(value);
    return [];
  };

  // Directories first, then alphabetical (case-insensitive).
  const sortFiles = (files: DG.FileInfo[]): void => {
    files.sort((a, b) => {
      if (a.isDirectory !== b.isDirectory)
        return a.isDirectory ? -1 : 1;
      return a.fileName.toLowerCase().localeCompare(b.fileName.toLowerCase());
    });
  };

  const renderFileLabel = (f: DG.FileInfo): HTMLElement => {
    const oh = DG.ObjectHandler.forEntity(f);
    return ui.divH([...(oh ? [oh.renderIcon(f.dart)] : []), ui.divText(f.fileName)], {style: {alignItems: 'center', gap: '4px'}});
  };

  // Lazily populates a group's children once. The promise is cached per Dart node handle so
  // that a manual expand and a programmatic restore never double-load the same group.
  const loadCache = new Map<unknown, Promise<void>>();
  const ensureLoaded = (group: DG.TreeViewGroup): Promise<void> => {
    let p = loadCache.get(group.dart);
    if (!p) {
      p = loadChildren(group);
      loadCache.set(group.dart, p);
    }
    return p;
  };

  const loadChildren = async (group: DG.TreeViewGroup): Promise<void> => {
    const removeLoader = addLoader(group);
    try {
      const files = await listContent(group.value);
      sortFiles(files);
      for (const f of files) {
        if (f.isDirectory) {
          const sub = group.group(renderFileLabel(f), f, false); // folders start collapsed
          stampFolder(sub.root, f.fileName);
        } else {
          const item = group.item(renderFileLabel(f), f);
          stampFile(item.root, f);
          ui.makeDraggable(item.root, {
            getDragObject: () => f,
            getDragCaption: () => f.friendlyName,
          });
        }
      }
    } catch (e) {
      grok.shell.error('Failed to list files. Check console for details');
      console.error(e);
    } finally {
      removeLoader();
    }
  };

  // Lazy loading: fires for any descendant group being expanded, at any depth.
  tree.onChildNodeExpanding.subscribe((n) => {
    ensureLoaded(n);
  });

  // Selection / activation events fire only for actual files.
  tree.onSelectedNodeChanged.subscribe((n) => {
    const v = n?.value;
    if (v instanceof DG.FileInfo && v.isFile)
      onSelected(v);
  });
  // onNodeEnter is the equivalent of double clicking a node.
  tree.onNodeEnter.subscribe((n) => {
    const v = n?.value;
    if (v instanceof DG.FileInfo && v.isFile)
      onDblClick(v);
  });

  // --- Expanded-state persistence ---------------------------------------------------------

  const collectExpanded = (group: DG.TreeViewGroup): ExpandedState => {
    const state: ExpandedState = {};
    for (const child of group.children) {
      if (child instanceof DG.TreeViewGroup && child.expanded) {
        const key = nodeKey(child);
        if (key != null)
          state[key] = collectExpanded(child);
      }
    }
    return state;
  };

  const saveExpandedState = (): void => {
    if (!localStorageKey)
      return;
    try {
      localStorage.setItem(localStorageKey, JSON.stringify(collectExpanded(tree)));
    } catch (e) {
      console.error(e);
    }
  };

  // Expands the groups described by `state`, in order, loading each level before descending.
  // A group that no longer exists stops expansion of everything downstream of it.
  const expandFromState = async (group: DG.TreeViewGroup, state: ExpandedState): Promise<void> => {
    for (const key of Object.keys(state)) {
      const child = group.children.find(
        (c): c is DG.TreeViewGroup => c instanceof DG.TreeViewGroup && nodeKey(c) === key);
      if (!child)
        continue; // group is gone — do not descend into this branch
      await ensureLoaded(child); // make sure nested groups exist before recursing
      child.expanded = true;
      await expandFromState(child, state[key]);
    }
  };

  const restoreExpandedState = async (): Promise<void> => {
    if (!localStorageKey)
      return;
    let saved: ExpandedState;
    try {
      const raw = localStorage.getItem(localStorageKey);
      if (!raw)
        return;
      saved = JSON.parse(raw);
    } catch (e) {
      console.error(e);
      return;
    }
    await expandFromState(tree, saved);
  };

  if (localStorageKey)
    tree.onChildNodeExpandedChanged.subscribe(() => saveExpandedState());

  // --- Initial population ------------------------------------------------------------------

  const processTree = async () => {
    const removeMainTreeLoader = addLoader(tree);
    try {
      const dataConnections = (await grok.dapi.connections.list())
        .filter((a) => a.dataSource == DG.DataSourceType.Files || DG.DataSourceType.fileDataSources.includes(a.dataSource))
        .map((ds) => ({con: ds, name: ds.friendlyName ?? ds.name}));

      removeMainTreeLoader();

      for (const dc of dataConnections) {
        const oh = DG.ObjectHandler.forEntity(dc.con);
        const g = tree.group(ui.divH([...(oh ? [oh.renderIcon(dc.con.dart)] : []), ui.divText(dc.name)], {style: {alignItems: 'center', gap: '4px'}}), dc.con, false);
        stampConnection(g.root, dc.name);
      }

      await restoreExpandedState();
    } catch (e) {
      removeMainTreeLoader();
      grok.shell.error('Something went wrong. Check console for details');
      console.error(e);
    }
  };

  processTree();
  tree.root.style.maxHeight = '300px';
  tree.root.style.minWidth = '200px';
  tree.root.style.minHeight = '50px';
  tree.root.classList.add('d4-tree-view-sticky');
  return tree;
}
