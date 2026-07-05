/** Hierarchical space chooser for dialogs: browse root spaces, drill into
 *  subspaces of any depth (loaded lazily on expand), create a subspace under
 *  the selected space (or a new root space), and select the target.
 *
 *  Selection is optional by design — the host dialog decides what "no space"
 *  means (Save As treats it as "save as a plain script in my namespace"). */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class SpacePicker {
  readonly root: HTMLElement;
  private readonly tree = DG.TreeViewGroup.tree();
  /** Space ids whose subspaces are already in the tree (lazy-load guard). */
  private readonly loaded = new Set<string>();
  private selection: DG.Project | null = null;

  /** Fires whenever the selected space changes (null = nothing selected). */
  onChanged: ((space: DG.Project | null) => void) | null = null;

  private constructor() {
    this.tree.root.classList.add('funcflow-space-picker-tree');
    this.tree.root.style.maxHeight = '220px';
    this.tree.root.style.minHeight = '50px';
    this.tree.root.style.width = '250px';
    this.tree.root.style.paddingBottom = '5px';
    this.tree.root.style.overflowY = 'auto';
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      this.select((node?.value as DG.Project) ?? null);
    });
    const newSpaceBtn = ui.button('New subspace…', () => void this.createSpaceDialog());
    ui.tooltip.bind(newSpaceBtn,
      'Create a subspace under the selected space, or a new root space when nothing is selected');
    this.root = ui.divV([this.tree.root, newSpaceBtn], 'funcflow-space-picker');
  }

  /** Builds a picker with the root spaces loaded; subspaces load on expand. */
  static async create(): Promise<SpacePicker> {
    const picker = new SpacePicker();
    let roots: DG.Project[] = [];
    try {
      roots = await grok.dapi.spaces.list();
    } catch { /* spaces may be unavailable on this server */ }
    for (const space of roots.filter((s) => s.friendlyName?.toLowerCase() !== 'admin'))
      picker.addSpaceNode(picker.tree, space);
    if (roots.length === 0)
      picker.tree.root.appendChild(ui.divText('No spaces are available on this server'));
    return picker;
  }

  get selected(): DG.Project | null { return this.selection; }

  private select(space: DG.Project | null): void {
    this.selection = space;
    this.onChanged?.(space);
  }

  private addSpaceNode(parent: DG.TreeViewGroup, space: DG.Project): DG.TreeViewGroup {
    const node = parent.group(space.friendlyName || space.name, space, false);
    ui.tooltip.bind(node.captionLabel, 'Click to select; expand to browse subspaces');
    node.onNodeExpanding.subscribe(() => void this.loadSubspaces(node, space));
    return node;
  }

  private async loadSubspaces(node: DG.TreeViewGroup, space: DG.Project): Promise<void> {
    if (this.loaded.has(space.id)) return;
    this.loaded.add(space.id);
    try {
      const children = await grok.dapi.spaces.id(space.id).children.filter('Project', false).list();
      for (const sub of children as DG.Project[])
        if (sub.isSpace)
          this.addSpaceNode(node, sub);
    } catch (e: any) {
      this.loaded.delete(space.id); // retry on the next expand
      grok.shell.error(`Could not load subspaces of "${space.friendlyName}": ${e?.message ?? e}`);
    }
  }

  private async createSpaceDialog(): Promise<void> {
    const parentNode = this.tree.currentItem instanceof DG.TreeViewGroup ? this.tree.currentItem : null;
    const parentSpace = (parentNode?.value as DG.Project) ?? null;
    const nameInput = ui.input.string('Name', {value: '',
      tooltipText: 'Name of the new space'});
    const where = parentSpace ? `Created under "${parentSpace.friendlyName}"` : 'Created as a new root space';

    ui.dialog({title: 'New Space'})
      .add(ui.divText(where))
      .add(nameInput)
      .onOK(async () => {
        const name = nameInput.value.trim();
        if (name === '') {
          grok.shell.warning('Give the space a name first');
          return;
        }
        try {
          let created: DG.Project;
          if (parentSpace != null) {
            if (await grok.dapi.spaces.id(parentSpace.id).subspaceExists(name)) {
              grok.shell.warning(`Subspace "${name}" already exists under "${parentSpace.friendlyName}"`);
              return;
            }
            created = await grok.dapi.spaces.id(parentSpace.id).addSubspace(name);
          } else {
            if (await grok.dapi.spaces.rootSpaceExists(name)) {
              grok.shell.warning(`Space "${name}" already exists`);
              return;
            }
            created = await grok.dapi.spaces.createRootSpace(name);
          }
          if (parentNode != null && parentSpace != null) {
            // Only add manually when the parent's children are already loaded —
            // otherwise the lazy load would bring the new subspace in twice.
            if (this.loaded.has(parentSpace.id))
              this.tree.currentItem = this.addSpaceNode(parentNode, created);
            parentNode.expanded = true;
          } else
            this.tree.currentItem = this.addSpaceNode(this.tree, created);
          this.select(created);
          grok.shell.info(`Space "${name}" created`);
        } catch (e: any) {
          // Server-side checks: CreateSpace privilege, name uniqueness, EDIT permission.
          grok.shell.error(`Could not create space: ${e?.message ?? e}`);
        }
      })
      .show();
  }
}
