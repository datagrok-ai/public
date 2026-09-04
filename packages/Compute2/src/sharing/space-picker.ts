import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Space chooser over a lazily-loaded tree (adapted from Flow's SpacePicker, chooser mode only). */
export class SpacePicker {
  readonly root: HTMLElement;
  readonly tree = DG.TreeViewGroup.tree();
  private readonly loaded = new Set<string>();
  private selection: DG.Project | null = null;

  onChanged: ((space: DG.Project | null) => void) | null = null;

  private constructor() {
    this.tree.root.style.maxHeight = '220px';
    this.tree.root.style.minHeight = '50px';
    this.tree.root.style.width = '100%';
    this.tree.root.style.overflowY = 'auto';
    this.tree.root.style.paddingBottom = '5px';
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      const v = node?.value;
      this.selection = v instanceof DG.Project && v.isSpace ? v : null;
      this.onChanged?.(this.selection);
    });
    this.root = ui.divV([this.tree.root]);
  }

  /** Builds a picker with the root spaces loaded; subspaces load on expand. */
  static async create(): Promise<SpacePicker> {
    const picker = new SpacePicker();
    let roots: DG.Project[] = [];
    try {
      roots = await grok.dapi.spaces.list();
    } catch {
      // spaces may be unavailable on this server
    }
    for (const space of roots.filter((s) => s.friendlyName?.toLowerCase() !== 'admin'))
      picker.addSpaceNode(picker.tree, space);
    if (roots.length === 0)
      picker.tree.root.appendChild(ui.divText('No spaces are available on this server'));
    return picker;
  }

  get selected(): DG.Project | null {
    return this.selection;
  }

  private addSpaceNode(parent: DG.TreeViewGroup, space: DG.Project): DG.TreeViewGroup {
    const name = space.friendlyName || space.name;
    const icon = document.createElement('i');
    icon.className = 'grok-icon fal fa-brackets-curly';
    const label = ui.divH([icon, ui.divText(name)], {style: {alignItems: 'center', gap: '4px'}});
    const node = parent.group(label, space, false);
    ui.tooltip.bind(node.captionLabel, 'Click to select; expand to browse subspaces');
    node.onNodeExpanding.subscribe(() => void this.loadChildren(node, space));
    return node;
  }

  private async loadChildren(node: DG.TreeViewGroup, space: DG.Project): Promise<void> {
    if (this.loaded.has(space.id)) return;
    this.loaded.add(space.id);
    try {
      const children = await grok.dapi.spaces.id(space.id).children.filter('Project', false).list();
      const seen = new Set<string>();
      const byName = (p: DG.Project) => (p.friendlyName || p.name).toLowerCase();
      const subs = children
        .filter((c): c is DG.Project => c instanceof DG.Project && c.isSpace)
        .filter((c) => {
          if (!c.id || seen.has(c.id)) return !c.id;
          seen.add(c.id);
          return true;
        })
        .sort((a, b) => byName(a).localeCompare(byName(b)));
      for (const sub of subs)
        this.addSpaceNode(node, sub);
    } catch (e: any) {
      this.loaded.delete(space.id); // retry on the next expand
      grok.shell.error(`Could not load subspaces of "${space.friendlyName}": ${e?.message ?? e}`);
    }
  }
}
