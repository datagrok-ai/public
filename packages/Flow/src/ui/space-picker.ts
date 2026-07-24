/** Hierarchical space browser used in two modes:
 *  - **Chooser** (default — the Save dialog): browse root spaces, drill into
 *    subspaces of any depth (loaded lazily on expand), create a subspace under
 *    the selected space (or a new root space), and select the target.
 *    Selection is optional by design — the host dialog decides what "no space"
 *    means (Save As treats it as "save as a plain script in my namespace").
 *  - **Content browser** (`showContent: true` — the toolbox Spaces tab): also
 *    lists each space's CONTENT — entities (flows, scripts, queries, …) and
 *    files stored in the space — as activatable/draggable rows. The
 *    "New subspace…" button is hidden (it's a browsing surface, not a
 *    save-target chooser). */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {setTid} from '../utils/test-ids';

export interface SpacePickerOptions {
  /** List each space's content (entities + stored files) alongside its
   *  subspaces, and hide the "New subspace…" button. Default: spaces only. */
  showContent?: boolean;
  /** Content mode: fired when a content row is double-clicked (or activated
   *  with Enter). Never fires for spaces — expanding is their interaction. */
  onEntityActivated?: (entity: DG.Entity) => void;
}

/** Display name for a content row (files show their full file name). */
export function spaceEntityName(e: DG.Entity): string {
  if (e instanceof DG.FileInfo) return e.fileName;
  return e.friendlyName || e.name || '';
}

export class SpacePicker {
  readonly root: HTMLElement;
  /** Exposed for tests (expanding a group triggers the lazy child load). */
  readonly tree = DG.TreeViewGroup.tree();
  /** Space ids whose children are already in the tree (lazy-load guard). */
  private readonly loaded = new Set<string>();
  private selection: DG.Project | null = null;
  private readonly options: SpacePickerOptions;

  /** Fires whenever the selected space changes (null = nothing selected). */
  onChanged: ((space: DG.Project | null) => void) | null = null;

  private constructor(options: SpacePickerOptions) {
    this.options = options;
    this.tree.root.classList.add('funcflow-space-picker-tree');
    if (options.showContent) {
      // The hosting pane (toolbox tab) owns width and scrolling.
      this.tree.root.style.width = '100%';
    } else {
      this.tree.root.style.maxHeight = '220px';
      this.tree.root.style.minHeight = '50px';
      this.tree.root.style.width = '250px';
      this.tree.root.style.overflowY = 'auto';
    }
    this.tree.root.style.paddingBottom = '5px';
    this.tree.onSelectedNodeChanged.subscribe((node) => {
      const v = node?.value;
      this.select(v instanceof DG.Project && v.isSpace ? v : null);
    });
    if (options.showContent) {
      // Double-click / Enter on a content row → hand the entity to the host.
      this.tree.onNodeEnter.subscribe((node) => {
        const v = node?.value;
        if (v instanceof DG.Project || !(v instanceof DG.Entity)) return;
        options.onEntityActivated?.(v);
      });
      this.root = ui.divV([this.tree.root], 'funcflow-space-picker');
    } else {
      const newSpaceBtn = ui.button('New subspace…', () => void this.createSpaceDialog());
      ui.tooltip.bind(newSpaceBtn,
        'Create a subspace under the selected space, or a new root space when nothing is selected');
      this.root = ui.divV([this.tree.root, newSpaceBtn], 'funcflow-space-picker');
    }
  }

  /** Builds a picker with the root spaces loaded; children load on expand. */
  static async create(options: SpacePickerOptions = {}): Promise<SpacePicker> {
    const picker = new SpacePicker(options);
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
    const name = space.friendlyName || space.name;
    // Space rows carry the Browse tree's space glyph — every other row type has
    // an icon, and it teaches what the toolbox tab's brackets-curly means.
    const icon = document.createElement('i');
    icon.className = 'grok-icon fal fa-brackets-curly funcflow-space-icon';
    const label = ui.divH([icon, ui.divText(name)], {style: {alignItems: 'center', gap: '4px'}});
    const node = parent.group(label, space, false);
    setTid(node.root, 'space', name);
    ui.tooltip.bind(node.captionLabel, this.options.showContent ?
      'Expand to browse subspaces and content' : 'Click to select; expand to browse subspaces');
    node.onNodeExpanding.subscribe(() => void this.loadChildren(node, space));
    return node;
  }

  /** One activatable/draggable row for a space's content entity (content mode).
   *  Files and functions drag onto the canvas exactly like toolbox rows — the
   *  canvas droppable already accepts DG.FileInfo and DG.Func drag objects. */
  private addContentNode(parent: DG.TreeViewGroup, entity: DG.Entity): void {
    const name = spaceEntityName(entity);
    const oh = DG.ObjectHandler.forEntity(entity);
    const label = ui.divH([...(oh ? [oh.renderIcon(entity.dart)] : []), ui.divText(name)],
      {style: {alignItems: 'center', gap: '4px'}});
    const item = parent.item(label, entity);
    setTid(item.root, 'space-item', name);
    const actionable = entity instanceof DG.Func || (entity instanceof DG.FileInfo && entity.isFile);
    if (actionable) {
      ui.tooltip.bind(item.root, `${name} — double-click or drag onto the canvas to use it`);
      ui.makeDraggable(item.root, {
        getDragObject: () => entity,
        getDragCaption: () => name,
      });
    }
  }

  private async loadChildren(node: DG.TreeViewGroup, space: DG.Project): Promise<void> {
    if (this.loaded.has(space.id)) return;
    this.loaded.add(space.id);
    try {
      const client = grok.dapi.spaces.id(space.id).children;
      // Empty types = every entity kind + the space's stored files; linked
      // entities included so a flow shared into the space shows up too.
      const children = this.options.showContent ?
        await client.filter('', true).list() :
        await client.filter('Project', false).list();
      const seen = new Set<string>();
      const unique = children.filter((c) => {
        if (!c.id || seen.has(c.id)) return !c.id;
        seen.add(c.id);
        return true;
      });
      const byName = (a: DG.Entity, b: DG.Entity): number =>
        spaceEntityName(a).toLowerCase().localeCompare(spaceEntityName(b).toLowerCase());
      const spaces = unique
        .filter((c): c is DG.Project => c instanceof DG.Project && c.isSpace).sort(byName);
      for (const sub of spaces)
        this.addSpaceNode(node, sub);
      if (this.options.showContent) {
        // Subspaces first, then entities, files last — mirrors the Files tree's
        // dirs-before-files reading order. Connections are skipped: every space
        // carries its own "<name> Files" storage connection (pure plumbing —
        // the stored files list as rows of their own), and a connection row
        // isn't actionable on the canvas anyway (its queries are entities).
        const rest = unique.filter((c) =>
          !(c instanceof DG.Project && c.isSpace) && !(c instanceof DG.DataConnection));
        const entities = rest.filter((c) => !(c instanceof DG.FileInfo)).sort(byName);
        const files = rest.filter((c) => c instanceof DG.FileInfo).sort(byName);
        for (const e of [...entities, ...files])
          this.addContentNode(node, e);
      }
    } catch (e: any) {
      this.loaded.delete(space.id); // retry on the next expand
      grok.shell.error(`Could not load content of "${space.friendlyName}": ${e?.message ?? e}`);
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
