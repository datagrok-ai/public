import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PieChartSettings, Sector, Subsector} from '../sparklines/piechart';
import {
  CONSTANTS,
  DEFAULTS,
  TAGS,
  defaultGeneralProps,
  defaultGroupProps,
  generalProps,
  groupProps,
  SectorType,
  subGroupProps,
  VlaaivisColumnMetadata
} from './constants';

import {MpoDesirabilityLineEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-line-editor';
import {NumericalDesirability} from '@datagrok-libraries/statistics/src/mpo/mpo';

class VlaaiVisManager {
  private settings: PieChartSettings;
  private gc: DG.GridColumn;
  private columns: DG.Column[];
  private tree: DG.TreeViewGroup;
  private metadataMap: Map<string, VlaaivisColumnMetadata> = new Map();

  constructor(settings: PieChartSettings, gc: DG.GridColumn) {
    this.settings = settings;
    this.gc = gc;
    this.columns = grok.shell.tv.dataFrame.columns.byNames(settings.columnNames);
    this.tree = ui.tree();
    this.tree.root.style.overflow = 'hidden';

    if (!this.settings.sectors || this.metadataMap.size === 0) {
      this.settings.sectors = this.initializeSectors(settings.columnNames, grok.shell.tv.dataFrame);
      this.gc.grid.invalidate();
    }
  }

  private initializeSectors(columnNames: string[], dataFrame: DG.DataFrame): {
    lowerBound: number;
    upperBound: number;
    sectors: Sector[];
    values: string | null;
  } {
    const groupMap = new Map<string, VlaaivisColumnMetadata[]>();
    const columns = dataFrame.columns.byNames(columnNames);

    for (const col of columns) {
      const tag = col.getTag(TAGS.VLAAIVIS_METADATA);
      if (!tag) continue;

      const meta: VlaaivisColumnMetadata = {...JSON.parse(tag), name: col.name};
      this.metadataMap.set(col.name, meta);

      if (!meta.groupName) continue;

      let group = groupMap.get(meta.groupName);
      if (!group) {
        group = [];
        groupMap.set(meta.groupName, group);
      }
      group[group.length] = meta;
    }

    const sectors: Sector[] = Array.from(groupMap, ([groupName, metas]) => ({
      name: groupName,
      sectorColor: metas[0].sectorColor ?? defaultGroupProps[CONSTANTS.SECTOR_COLOR_PROPERTY],
      subsectors: metas.map(({name, weight, line, min, max}) => ({name, functionType: 'numerical' as const, weight, line, min, max}))
    }));

    return {
      lowerBound: defaultGeneralProps[CONSTANTS.LOWER_BOUND],
      upperBound: defaultGeneralProps[CONSTANTS.UPPER_BOUND],
      sectors,
      values: null
    };
  }

  private findSectorOrSubsector(name: string): { entity: Sector | Subsector; type: SectorType } | null {
    const sectors = this.settings.sectors;
    if (!sectors) return null;

    for (const sector of sectors.sectors) {
      if (sector.name === name) {
        const result = {entity: sector, type: SectorType.SECTOR};
        return result;
      }

      const subsector = sector.subsectors.find((sub) => sub.name === name);
      if (subsector) {
        const result = {entity: subsector, type: SectorType.SUBSECTOR};
        return result;
      }
    }
    return null;
  }

  private getSectorProperty(name: string, propertyName: keyof (Sector | Subsector)) {
    const result = this.findSectorOrSubsector(name);
    if (result)
      return result.entity[propertyName];

    return name === '' ? (this.settings.sectors as any)[propertyName] : undefined;
  }

  private updateSectorProperty(name: string, propertyName: keyof (Sector | Subsector), value: any): void {
    const columnMap = new Map<string, DG.Column>(this.columns.map((col) => [col.name, col]));
    const sectorOrSubsector = this.findSectorOrSubsector(name);

    if (sectorOrSubsector) {
      sectorOrSubsector.entity[propertyName] = value;
      propertyName === CONSTANTS.SECTOR_COLOR_PROPERTY ?
        this.updateSectorColorTags(name, value) :
        this.updateColumnTag(name, propertyName, value, columnMap);
    } else if (name === '') {
      (this.settings.sectors as any)[propertyName] = value;
    }
  }

  private updateSectorColorTags(name: string, value: any): void {
    for (const [colName, meta] of this.metadataMap) {
      if (meta.groupName === name) {
        meta.sectorColor = value;
        const col = this.columns.find((c) => c.name === colName);
        if (col)
          col.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
      }
    }
  }

  private updateColumnTag(
    name: string, propertyName: keyof (Sector | Subsector), value: any, columnMap: Map<string, DG.Column>
  ): void {
    const column = columnMap.get(name);
    const meta = this.metadataMap.get(name) || this.metadataMap.set(name, {} as VlaaivisColumnMetadata).get(name);
    if (!meta || !column) return;
    meta[propertyName] = value;
    column.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
  }

  private createInputForProperty(property: any, name: string): HTMLElement {
    const prop = DG.Property.fromOptions(property.property);
    const input = DG.InputBase.forProperty(prop, {});
    input.enabled = property.property.enabled ?? true;
    const propName = property.property.name;

    const existingValue = this.getSectorProperty(name, propName);

    if (existingValue !== undefined) {
      input.value = existingValue;
    } else {
      const column = this.columns.find((c) => c.name === name);
      switch (propName) {
      case 'min':
        input.value = +column!.min.toFixed(1);
        break;
      case 'max':
        input.value = +column!.max.toFixed(1);
        break;
      default:
        input.value = property.object?.[propName];
        break;
      }
    }

    input.onChanged.subscribe((value) => {
      this.updateSectorProperty(name, propName, value);
      // TODO: Handle min/max change if needed
      this.gc.grid.invalidate();
    });

    return input.root;
  }

  private makeItemDraggable(item: DG.TreeViewNode<any>): void {
    item.root.onmouseenter = (e) => ui.tooltip.show('To visualise, drag item into existing group', e.x, e.y);
    item.root.onmouseleave = () => ui.tooltip.hide();

    ui.makeDraggable(item.root, {
      getDragObject: () => item,
      getDragCaption: () => `You are dragging ${item.text}`
    });
  }

  private makeGroupDroppable(groupNode: DG.TreeViewGroup): void {
    ui.makeDroppable(groupNode.root, {
      acceptDrop: () => true,
      doDrop: (draggedItem: any) => {
        const itemText = draggedItem.text;
        this.removeItemFromOriginalGroup(itemText);
        this.addItemToNewGroup(groupNode, itemText);
        this.gc.grid.invalidate();
      }
    });
  }

  private removeItemFromOriginalGroup(itemText: string): void {
    const group = this.tree.items.find((g) => g.text === itemText);
    if (group) {
      group.remove();
      const originalGroup = this.settings.sectors?.sectors.find((sector) => sector.name === group.parent.text);
      if (originalGroup)
        originalGroup.subsectors = originalGroup.subsectors.filter((subsector) => subsector.name !== itemText);
    }
  }

  private addItemToNewGroup(groupNode: DG.TreeViewGroup, itemText: string): void {
    const newGroup = this.getOrCreateNewGroup(groupNode.text);
    const column = this.columns.find((col) => col.name === itemText);

    this.updateColumnTags(column);

    if (newGroup && !newGroup.subsectors.some((subsector) => subsector.name === itemText))
      newGroup.subsectors = [...newGroup.subsectors, this.createSubsector(column)];

    if (!column) return;

    const {name} = column;
    const meta = this.metadataMap.get(name) || this.metadataMap.set(name, {} as VlaaivisColumnMetadata).get(name);
    const groupMeta = [...this.metadataMap.values()].find((m) => m.groupName === groupNode.text);
    if (!meta) return;

    meta.groupName = groupNode.text;
    meta.sectorColor = groupMeta?.sectorColor;
    column.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
    const newItem = groupNode.item(itemText);
    this.makeItemDraggable(newItem);
  }

  private getOrCreateNewGroup(groupName: string): Sector | null {
    if (groupName === '')
      return null;

    let newGroup = this.settings.sectors?.sectors.find((sector) => sector.name === groupName);
    if (!newGroup) {
      const defaultColor = groupProps[0].object['sectorColor'];
      const sectorColor = this.metadataMap.get(groupName ?? '')?.sectorColor ?? defaultColor;

      newGroup = {
        name: groupName,
        sectorColor: sectorColor,
        subsectors: []
      };
      this.settings.sectors!.sectors = [...this.settings.sectors!.sectors, newGroup];
    }
    return newGroup;
  }

  private createSubsector(column?: DG.Column): Subsector {
    const name = column?.name ?? '';
    const meta = name ? this.metadataMap.get(name) : {} as VlaaivisColumnMetadata;
    const metadataWeight = meta?.weight ?? undefined;
    const weight = typeof metadataWeight === 'number' ? metadataWeight : +this.generateRandomNumber().toFixed(1);

    return {
      name,
      functionType: 'numerical',
      weight,
      line: meta?.line ?? []
    };
  }

  private updateColumnTags(column: DG.Column | undefined): void {
    if (!column) return;

    const defaultWeight = this.generateRandomNumber().toFixed(1);
    const tag = column.getTag(TAGS.VLAAIVIS_METADATA);
    const meta: VlaaivisColumnMetadata = tag ? JSON.parse(tag) : {};

    meta.min ??= +column.min.toFixed(1);
    meta.max ??= +column.max.toFixed(1);
    meta.weight ??= +defaultWeight;

    column.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
    this.metadataMap.set(column.name, meta);
  }

  private initializeTreeGroups(): Map<string, DG.Column[]> {
    const groupMap = new Map<string, DG.Column[]>();

    this.columns.forEach((col) => {
      const metadata = this.metadataMap.get(col.name);
      const groupName = metadata?.groupName;
      if (groupName) {
        const existingColumns = groupMap.get(groupName) ?? [];
        groupMap.set(groupName, [...existingColumns, col]);
      }
    });

    return groupMap;
  }

  private createLineEditor(nodeText: string, container: HTMLElement): void {
    const {line = [], min, max, weight} = this.metadataMap.get(nodeText) ?? {};
    const column = this.columns.find((c) => c.name === nodeText)!;

    const lineProp: NumericalDesirability = {
      functionType: 'numerical',
      line,
      min: min ?? +column.min.toFixed(1),
      max: max ?? +column.max.toFixed(1),
      weight: weight ?? +DEFAULTS.WEIGHT,
    };

    const lineEditor = new MpoDesirabilityLineEditor(lineProp, 200, 80);
    lineEditor.root.classList.add('mpo-line-editor');
    this.updateSectorProperty(nodeText, 'line' as keyof (Sector | Subsector), lineEditor.line);

    lineEditor.onChanged.subscribe(() => {
      this.updateSectorProperty(nodeText, 'line' as keyof (Sector | Subsector), lineEditor.line);
      this.gc.grid.invalidate();
    });

    const oldEditor = container.querySelector('.mpo-line-editor');
    oldEditor?.remove();

    container.appendChild(lineEditor.root);
  }

  public createTreeGroup(): HTMLElement {
    const inputs = ui.divV([], {style: {alignItems: 'end'}});
    const groupMap = this.initializeTreeGroups();
    const untaggedColumns = this.getUntaggedColumns();

    groupMap.forEach((groupColumns, groupName) => {
      const groupNode = this.tree.group(groupName);
      this.configureGroupNode(groupNode);

      groupColumns.forEach((col) => {
        const colNode = groupNode.item(col.name);
        this.makeItemDraggable(colNode);
      });
    });

    this.tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
      const nodeText = node?.text;
      const isValidNode = node.parent.text !== '';
      const column = this.columns.find((c) => c.name === nodeText);

      if (!isValidNode || !column) return;

      this.updateInputs(inputs, nodeText);
      this.createLineEditor(nodeText, inputs);
    });


    this.tree.onNodeContextMenu.subscribe((args: any) => {
      const menu: DG.Menu = args.args.menu;
      const node: DG.TreeViewNode = args.args.item;
      if (node instanceof DG.TreeViewGroup) {
        menu.item('Edit...', () => this.editGroup(node));
        menu.item('Delete', () => this.deleteGroup(node));
      } else { menu.item('Add...', () => this.addGroup(node)); }
    });

    this.makeUntaggedColumnsDraggable(untaggedColumns);

    const generalInputs = ui.divV(generalProps.map((prop) => this.createInputForProperty(prop, '')));
    const resultingDiv = this.createResultingDiv(inputs, generalInputs);
    return resultingDiv;
  }

  private configureGroupNode(groupNode: DG.TreeViewGroup): void {
    const colorPicker = this.createColorPicker(groupNode.text);
    const referenceNode = groupNode.root.getElementsByClassName('d4-tree-view-group-label')[0];
    referenceNode.insertAdjacentElement('beforebegin', colorPicker);

    groupNode.expanded = true;
    this.makeGroupDroppable(groupNode);
  }

  private createColorPicker(groupName: string): HTMLElement {
    const colorPicker = this.createInputForProperty(groupProps[0], groupName)
      .getElementsByClassName('ui-input-options')[0] as HTMLElement;
    (colorPicker as HTMLElement).style.cssText += ('margin: 5px!important');
    return colorPicker;
  }

  private updateInputs(inputs: HTMLElement, nodeText: string): void {
    ui.empty(inputs);
    if (this.settings.columnNames.includes(nodeText)) {
      subGroupProps.forEach((prop) => {
        inputs.appendChild(this.createInputForProperty(prop, nodeText));
      });
    }
  }

  private createResultingDiv(inputs: HTMLDivElement, generalInputs: HTMLDivElement): HTMLElement {
    generalInputs.style.marginRight = '10px';
    const tree = ui.divH([this.tree.root, inputs], 'ui-form');
    const container = ui.divV([generalInputs, tree]);
    container.style.marginTop = '10px';
    return container;
  }

  private addGroup(node: DG.TreeViewNode) {
    const nameInput = ui.input.string('Name');
    ui.dialog('New group')
      .add(nameInput)
      .onOK(() => {
        const newGroup = this.tree.group(nameInput.value, null, false, this.findLastGroupIndex());
        this.configureGroupNode(newGroup);
        this.removeItemFromOriginalGroup(node.text);
        this.addItemToNewGroup(newGroup, node.text);
        this.makeItemDraggable(node);
        this.gc.grid.invalidate();
      })
      .show();
  }

  private editGroup(group: DG.TreeViewGroup) {
    const nameInput = ui.input.string('Name');
    ui.dialog('Edit')
      .add(nameInput)
      .onOK(() => {
        this.columns = this.columns.map((col) => {
          const meta = this.metadataMap.get(col.name);
          const groupName = meta?.groupName;
          if (groupName === group.text && meta) {
            meta.groupName = nameInput.stringValue;
            col.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
          }
          return col;
        });

        const sector = this.settings.sectors?.sectors.find((s) => s.name === group.text);
        if (sector)
          sector.name = nameInput.stringValue;

        group.text = nameInput.stringValue;
      })
      .show();
  }

  private deleteGroup(group: DG.TreeViewGroup) {
    const groupChildren = group.children;
    group.remove();

    if (!groupChildren) return;

    const sectors = this.settings.sectors!.sectors;
    const groupName = group.text;

    for (let i = 0; i < groupChildren.length; i++) {
      const child = groupChildren[i];
      if (!(child instanceof DG.TreeViewNode)) continue;

      const columnName = child.text;
      const column = this.columns.find((col) => col.name === columnName);
      if (!column) continue;

      const meta = this.metadataMap.get(columnName);
      if (!meta) continue;

      delete meta.groupName;
      delete meta.sectorColor;

      column.setTag(TAGS.VLAAIVIS_METADATA, JSON.stringify(meta));
      this.makeItemDraggable(this.tree.item(columnName));
    }

    this.settings.sectors!.sectors = sectors.filter((sector) => sector.name !== groupName);
    this.gc.grid.invalidate();
  }

  private getUntaggedColumns(): DG.Column[] {
    return this.columns.filter((col) => !(this.metadataMap.get(col.name)?.groupName));
  }

  private makeUntaggedColumnsDraggable(columns: DG.Column[]): void {
    columns.forEach((col) => {
      const colNode = this.tree.item(col.name);
      this.makeItemDraggable(colNode);
    });
  }

  private generateRandomNumber(): number {
    return Math.random();
  }

  private findLastGroupIndex(): number {
    return this.tree.items.reduce((count, item) =>
      count + (item instanceof DG.TreeViewGroup ? 1 : 0), 0);
  }
}

export {VlaaiVisManager};
