/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TreeViewGroup} from 'datagrok-api/dg';
import {BitwiseOp, ScaffoldTreeViewer, ensureNodeBitset, getMol, isOrphans, value} from '../widgets/scaffold-tree';

const _funcs = new Map<string, DG.Func>();

function aiFunc(name: string, signature: string, description: string,
  run: (viewer: ScaffoldTreeViewer, ...args: any[]) => Promise<any>): DG.Func {
  let f = _funcs.get(name);
  if (!f) {
    f = grok.functions.register({signature, isAsync: true, namespace: 'ChemScaffoldAI',
      tags: 'chem-scaffold-tree-ai-function', options: {description},
      run: async (widget: any, ...args: any[]) => run(viewerOf(widget), ...args)});
    _funcs.set(name, f);
  }
  return f;
}

const viewerByWrapper = new WeakMap<object, ScaffoldTreeViewer>();

export function attachScaffoldTreeAi(viewer: ScaffoldTreeViewer): void {
  const w = (viewer.dart ? DG.toJs(viewer.dart) : null) ?? viewer;
  viewerByWrapper.set(viewer, viewer);
  if (w !== viewer)
    viewerByWrapper.set(w, viewer);
}

function viewerOf(widget: any): ScaffoldTreeViewer {
  for (const candidate of [widget, widget?.jsView, widget?.jsViewer, widget?.dart ? DG.toJs(widget.dart) : null]) {
    const viewer = candidate != null ? viewerByWrapper.get(candidate) : null;
    if (viewer != null)
      return viewer;
  }
  for (const v of Array.from(grok.shell.tv?.viewers ?? [])) {
    const viewer = viewerByWrapper.get(v) ?? (v instanceof ScaffoldTreeViewer ? v : null);
    if (viewer != null)
      return viewer;
  }
  throw new Error('Target the Scaffold Tree viewer (pass its widget ref from list_view_widgets)');
}

type Scaffold = {group: TreeViewGroup, error?: undefined} | {group?: undefined, error: string};

function scaffoldAt(viewer: ScaffoldTreeViewer, path: string): Scaffold {
  let node: any = path ? viewer.tree : null;
  for (const index of (path ?? '').split('.'))
    node = node?.children?.[Number(index)];
  return node ? {group: node as TreeViewGroup} :
    {error: `No scaffold '${path}' - call getScaffoldTreeState for the current paths`};
}

function structureAt(viewer: ScaffoldTreeViewer, path: string): Scaffold {
  const found = scaffoldAt(viewer, path);
  return found.group && isOrphans(found.group) ? {error: `'${path}' is an orphans folder, not a scaffold`} : found;
}

async function matchedScaffoldAt(viewer: ScaffoldTreeViewer, path: string): Promise<Scaffold> {
  const found = scaffoldAt(viewer, path);
  if (found.group && !isOrphans(found.group))
    await ensureNodeBitset(viewer, found.group, false);
  return found;
}

function filterState(viewer: ScaffoldTreeViewer) {
  return {
    checkedCount: viewer.checkedNodes.length,
    filterOperation: viewer.bitOperation,
    filteredRows: viewer.dataFrame?.filter.trueCount ?? null,
  };
}

function moleculeError(molStr: string): string | null {
  if (!molStr)
    return 'Pass smiles - a SMILES or molblock';
  const mol = getMol(molStr);
  if (mol === null)
    return `Cannot parse '${molStr}' as a molecule`;
  mol.delete();
  return null;
}

export function scaffoldTreeBriefing(viewer: ScaffoldTreeViewer): string {
  return `Scaffold tree over ${viewer.moleculeColumnName || 'no molecule column'}: ` +
    `${viewer.tree.items.length} scaffolds, ${viewer.checkedNodes.length} of them filtering with ${viewer.bitOperation}. ` +
    'Act on it through this widget\'s functions (list_view_functions with "scaffold", targeting this widget ' +
    'ref) - never through datagrok-exec, they take a widget argument only the tool call supplies. ' +
    'getScaffoldTreeState first (it gives the scaffold paths the others take), then addScaffold, editScaffold, ' +
    'deleteScaffold, setScaffoldColor, filterByScaffold, setFilterOperation, clearScaffoldFilters, ' +
    'generateScaffoldTree, loadScaffoldTree, saveScaffoldTree, clearScaffoldTree, selectRowsForScaffold.';
}

export function scaffoldTreeFunctions(): DG.Func[] {
  return [
    aiFunc('getScaffoldTreeState', 'dynamic getScaffoldTreeState(widget widget)',
      'The tree as it is saved: every scaffold with its structure, checked (filtering) and isNot (excluding) flags, colors and label, nested through child_nodes. Call this first - the other functions address a scaffold by its path here: the index of the root scaffold, then the indices through child_nodes, e.g. 0 or 0.2',
      async (viewer) => {
        const blocked = viewer.generateBlockedReason();
        return {
          moleculeColumn: viewer.moleculeColumnName || null,
          scaffoldCount: viewer.tree.items.length,
          ...filterState(viewer),
          totalRows: viewer.dataFrame?.rowCount ?? null,
          ...(blocked ? {generateBlocked: blocked} : {}),
          scaffolds: JSON.parse(viewer.treeEncode, (key, node) => key === 'orphansBitset' ? undefined : node)
            .filter((node: any) => typeof node !== 'string'),
        };
      }),

    aiFunc('addScaffold', 'dynamic addScaffold(widget widget, string smiles, string parent)',
      'Add a scaffold - a SMILES or a molblock. parent is optional: the path to nest the new scaffold under, omitted for a root scaffold. A nested scaffold must contain its parent as a substructure',
      async (viewer, smiles: string, parent?: string) => {
        const invalid = moleculeError(smiles);
        if (invalid)
          return {success: false, error: invalid};

        let parentGroup = viewer.tree;
        if (parent) {
          const {group, error} = structureAt(viewer, parent);
          if (!group)
            return {success: false, error};
          const parentSmiles = value(group).smiles;
          if (!ScaffoldTreeViewer.validateNodes(smiles, parentSmiles))
            return {success: false, error: `'${smiles}' does not contain '${parentSmiles}' - a nested scaffold must be a superstructure of its parent`};
          parentGroup = group;
        }

        const child = await viewer.addScaffoldNode(smiles, parentGroup);
        if (child === null)
          return {success: false, error: 'No molecule column - the scaffold tree has nothing to match against'};

        const index = parentGroup.children.indexOf(child);
        return {success: true, path: parent ? `${parent}.${index}` : `${index}`};
      }),

    aiFunc('editScaffold', 'dynamic editScaffold(widget widget, string scaffold, string smiles)',
      'Replace the structure of the scaffold at a path. The new structure must still contain the parent scaffold and be contained in every child scaffold',
      async (viewer, scaffold: string, smiles: string) => {
        const invalid = moleculeError(smiles);
        if (invalid)
          return {success: false, error: invalid};
        const {group, error} = structureAt(viewer, scaffold);
        if (!group)
          return {success: false, error};
        const conflict = viewer.validateScaffoldEdit(group, smiles);
        if (conflict)
          return {success: false, error: conflict};
        await viewer.applyScaffoldEdit(group, smiles);
        return {success: true};
      }),

    aiFunc('deleteScaffold', 'dynamic deleteScaffold(widget widget, list scaffolds)',
      'Delete scaffolds by path, together with everything nested under them. This cannot be undone',
      async (viewer, scaffolds: string[]) => {
        const found = scaffolds.map((path) => scaffoldAt(viewer, path));
        const missing = found.find((x) => x.error);
        if (missing)
          return {success: false, error: missing.error};
        for (const x of found)
          viewer.deleteScaffoldNode(x.group!);
        return {success: true, scaffoldCount: viewer.tree.items.length};
      }),

    aiFunc('setScaffoldColor', 'dynamic setScaffoldColor(widget widget, string scaffold, string color)',
      'Color the scaffold at a path and highlight it in the molecules that contain it; nested scaffolds inherit the color. color is a hex string like #e66465 - omit it to turn the color off',
      async (viewer, scaffold: string, color?: string) => {
        const {group, error} = structureAt(viewer, scaffold);
        if (!group)
          return {success: false, error};
        const html = color ? color.trim() : null;
        if (html !== null && !/^#[0-9a-f]{6}$/i.test(html))
          return {success: false, error: 'color must be a hex string like #e66465, or omitted to clear the color'};
        viewer.setScaffoldColor(group, html);
        return {success: true};
      }),

    aiFunc('filterByScaffold', 'dynamic filterByScaffold(widget widget, string scaffold, string mode)',
      'Filter the table by the scaffold at a path. mode is optional: only (default - filter by this scaffold alone), add (add it to the scaffolds already filtering), exclude (add it, keeping the molecules that do NOT contain it), remove (stop filtering by it). Several added scaffolds combine with the AND / OR set by setFilterOperation',
      async (viewer, scaffold: string, mode?: string) => {
        const m = (mode ?? 'only').toLowerCase();
        if (!['only', 'add', 'exclude', 'remove'].includes(m))
          return {success: false, error: 'mode must be one of: only, add, exclude, remove'};

        const {group, error} = await matchedScaffoldAt(viewer, scaffold);
        if (!group)
          return {success: false, error};
        if (m === 'only')
          viewer.makeNodeActiveAndFilter(group);
        else {
          group.checked = m !== 'remove';
          if (m !== 'remove')
            viewer.setNotBitOperation(group, m === 'exclude');
          viewer.updateFilters();
        }
        return {success: true, ...filterState(viewer)};
      }),

    aiFunc('setFilterOperation', 'dynamic setFilterOperation(widget widget, string operation)',
      'How several filtering scaffolds combine: OR keeps the molecules matching any of them, AND only those matching all of them',
      async (viewer, operation: string) => {
        const op = (operation ?? '').toUpperCase();
        if (op !== BitwiseOp.AND && op !== BitwiseOp.OR)
          return {success: false, error: 'operation must be AND or OR'};
        viewer.getProperty('bitOperation')!.set(viewer, op);
        return {success: true, ...filterState(viewer)};
      }),

    aiFunc('clearScaffoldFilters', 'dynamic clearScaffoldFilters(widget widget)',
      'Stop filtering by every scaffold and bring all the rows back. The scaffolds themselves are kept - use clearScaffoldTree to delete them',
      async (viewer) => {
        viewer.clearFilters();
        return {success: true};
      }),

    aiFunc('generateScaffoldTree', 'dynamic generateScaffoldTree(widget widget)',
      'Build the scaffold tree from the molecule column, replacing whatever the tree holds now. Only works for small tables and may take a while',
      async (viewer) => {
        const blocked = viewer.generateBlockedReason();
        if (blocked)
          return {success: false, error: `${blocked}. Add scaffolds with addScaffold or load a saved tree instead`};
        await viewer.generateTree();
        return {success: true, scaffoldCount: viewer.tree.items.length};
      }),

    aiFunc('loadScaffoldTree', 'dynamic loadScaffoldTree(widget widget, string source)',
      'Load a saved scaffold tree, replacing the current one. source is either the path of a .tree file (e.g. System:AppData/Chem/demo.tree) or the tree JSON itself',
      async (viewer, source: string) => {
        const text = (source ?? '').trim();
        if (!text)
          return {success: false, error: 'Pass source - the path of a .tree file or the tree JSON'};

        const json = text.startsWith('[') ? text : await grok.dapi.files.readAsText(text);
        if (!Array.isArray(JSON.parse(json)))
          return {success: false, error: 'A scaffold tree is a JSON array of nodes'};

        await viewer.loadTreeStr(json);
        viewer.saveTreeEncode();
        return {success: true, scaffoldCount: viewer.tree.items.length};
      }),

    aiFunc('saveScaffoldTree', 'dynamic saveScaffoldTree(widget widget, string name)',
      'Save the tree as a .tree file in Downloads, the same way the Save tree command does. name is optional - the file name without the extension, e.g. my-tree. Load it back with loadScaffoldTree',
      async (viewer, name?: string) => {
        if (viewer.tree.children.length === 0)
          return {success: false, error: 'The tree is empty - nothing to save'};
        return {success: true, file: viewer.downloadTree(name)};
      }),

    aiFunc('clearScaffoldTree', 'dynamic clearScaffoldTree(widget widget)',
      'Delete every scaffold in the tree, drop the scaffold colors and bring all the rows back. This cannot be undone',
      async (viewer) => {
        viewer.deleteAllTrees();
        return {success: true};
      }),

    aiFunc('selectRowsForScaffold', 'dynamic selectRowsForScaffold(widget widget, string scaffold, bool select)',
      'Select the rows whose molecules contain the scaffold at a path, adding them to the current selection; pass select false to deselect them. Selection does not hide anything - use filterByScaffold to filter',
      async (viewer, scaffold: string, select?: boolean) => {
        const {group, error} = await matchedScaffoldAt(viewer, scaffold);
        if (!group)
          return {success: false, error};
        viewer.selectTableRows(group, select !== false);
        return {success: true, selectedRows: viewer.dataFrame?.selection.trueCount ?? null};
      }),
  ];
}
