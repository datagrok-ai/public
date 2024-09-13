import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import { PieChartSettings, Sector, Subsector } from '../sparklines/piechart';
import { CONSTANTS, defaultGeneralProps, defaultGroupProps, DEFAULTS, generalProps, groupProps, SectorType, subGroupProps, TAGS } from './constants';

class VlaaiVisManager {
  private settings: PieChartSettings;
  private gc: DG.GridColumn;
  private columns: DG.Column[];
  private tree: DG.TreeViewGroup;

  constructor(settings: PieChartSettings, gc: DG.GridColumn) {
    this.settings = settings;
    this.gc = gc;
    this.columns = grok.shell.tv.dataFrame.columns.byNames(settings.columnNames);
    this.tree = ui.tree();
    this.tree.root.style.overflow = 'hidden';

    if (!this.settings.sectors) {
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
    const groupMap = new Map<string, DG.Column[]>();
    const columns = dataFrame.columns.byNames(columnNames);

    columns.forEach(col => {
      const groupName = col.getTag(TAGS.GROUP_NAME);
      if (groupName) {
        const existingColumns = groupMap.get(groupName) ?? [];
        groupMap.set(groupName, [...existingColumns, col]);
        this.updateColumnTags(col);
      }
    });

    const sectors = Array.from(groupMap.entries()).map(([groupName, columns]) => ({
      name: groupName,
      sectorColor: columns[0].getTag(TAGS.SECTOR_COLOR) ?? defaultGroupProps[CONSTANTS.SECTOR_COLOR_PROPERTY],
      subsectors: columns.map(col => ({
        name: col.name,
        low: parseFloat(col.getTag(TAGS.LOW)) ?? DEFAULTS.LOW,
        high: parseFloat(col.getTag(TAGS.HIGH)) ?? DEFAULTS.HIGH,
        weight: parseFloat(col.getTag(TAGS.WEIGHT)) ?? parseFloat(this.generateRandomNumber().toFixed(1)),
      }))
    }));

    return {
      lowerBound: defaultGeneralProps[CONSTANTS.LOWER_BOUND],
      upperBound: defaultGeneralProps[CONSTANTS.UPPER_BOUND],
      sectors,
      values: null
    };
  }

  private findSectorOrSubsector(name: string): { entity: Sector | Subsector; type: SectorType } | null {
    const cache = new Map<string, { entity: Sector | Subsector; type: SectorType } | null>();

    if (cache.has(name)) return cache.get(name)!;

    const sectors = this.settings.sectors;
    if (!sectors) return null;

    for (const sector of sectors.sectors) {
      if (sector.name === name) {
        const result = { entity: sector, type: SectorType.SECTOR };
        cache.set(name, result);
        return result;
      }

      const subsector = sector.subsectors.find(sub => sub.name === name);
      if (subsector) {
        const result = { entity: subsector, type: SectorType.SUBSECTOR };
        cache.set(name, result);
        return result;
      }
    }

    cache.set(name, null);
    return null;
  }

  private getSectorProperty(name: string, propertyName: keyof (Sector | Subsector)) {
    const result = this.findSectorOrSubsector(name);
    return result ? result.entity[propertyName] : (name === '' ? (this.settings.sectors as any)[propertyName] : undefined);
  }

  private updateSectorProperty(name: string, propertyName: keyof (Sector | Subsector), value: any): void {
    const columnMap = new Map<string, DG.Column>(this.columns.map(col => [col.name, col]));
    const sectorOrSubsector = this.findSectorOrSubsector(name);

    if (sectorOrSubsector) {
      sectorOrSubsector.entity[propertyName] = value;
      propertyName === CONSTANTS.SECTOR_COLOR_PROPERTY 
        ? this.updateSectorColorTags(name, value)
        : this.updateColumnTag(name, propertyName, value, columnMap);
    } else if (name === '')
      (this.settings.sectors as any)[propertyName] = value;
  }

  private updateSectorColorTags(name: string, value: any): void {
    this.columns
      .filter(col => col.getTag(TAGS.GROUP_NAME) === name)
      .forEach(col => col.setTag(`${CONSTANTS.TAG_PREFIX}${CONSTANTS.SECTOR_COLOR_PROPERTY}`, value));
  }

  private updateColumnTag(name: string, propertyName: keyof (Sector | Subsector), value: any, columnMap: Map<string, DG.Column>): void {
    const column = columnMap.get(name);
    column?.setTag(`${CONSTANTS.TAG_PREFIX}${propertyName}`, value);
  }

  private createInputForProperty(property: any, name: string): HTMLElement {
    const prop = DG.Property.fromOptions(property.property);
    const input = DG.InputBase.forProperty(prop, {});
    const value = this.getSectorProperty(name, property.property.name as keyof (Sector | Subsector));
    input.value = value !== undefined ? value : property.object[property.property.name];
    
    input.onChanged.subscribe((value) => {
      this.updateSectorProperty(name, property.property.name as keyof (Sector | Subsector), value);
      this.gc.grid.invalidate();
    });

    return input.root;
  }

  private makeItemDraggable(item: DG.TreeViewNode<any>): void {
    item.root.onmouseenter = (e) => ui.tooltip.show('To visualise, drag item into existing group', e.x, e.y);
    item.root.onmouseleave = (e) => ui.tooltip.hide();

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
    const group = this.tree.items.find(g => g.text === itemText);
    if (group) {
      group.remove();
      const originalGroup = this.settings.sectors?.sectors.find(sector => sector.name === group.parent.text);
      if (originalGroup) {
        originalGroup.subsectors = originalGroup.subsectors.filter(subsector => subsector.name !== itemText);
      }
    }
  }

  private addItemToNewGroup(groupNode: DG.TreeViewGroup, itemText: string): void {
    const newGroup = this.getOrCreateNewGroup(groupNode.text);
    const column = this.columns.find(col => col.name === itemText);

    this.updateColumnTags(column);

    if (newGroup && !newGroup.subsectors.some(subsector => subsector.name === itemText))
      newGroup.subsectors = [...newGroup.subsectors, this.createSubsector(column)];

    if (column)
      column.setTag(TAGS.GROUP_NAME, groupNode.text);

    const newItem = groupNode.item(itemText);
    this.makeItemDraggable(newItem);
  }

  private getOrCreateNewGroup(groupName: string): Sector | null {
    if (groupName === '')
      return null;

    let newGroup = this.settings.sectors?.sectors.find(sector => sector.name === groupName);
    if (!newGroup) {
      newGroup = {
        name: groupName,
        sectorColor: this.columns.find(col => col.name === groupName)?.getTag(TAGS.SECTOR_COLOR) ?? groupProps[0].object['sectorColor'],
        subsectors: []
      };
      this.settings.sectors!.sectors = [...this.settings.sectors!.sectors, newGroup];
    }
    return newGroup;
  }

  private createSubsector(column: DG.Column | undefined): Subsector {
    return {
      name: column?.name ?? '',
      low: parseFloat(column?.getTag(TAGS.LOW) ?? DEFAULTS.LOW),
      high: parseFloat(column?.getTag(TAGS.HIGH) ?? DEFAULTS.HIGH),
      weight: parseFloat(column?.getTag(TAGS.WEIGHT) ?? this.generateRandomNumber().toFixed(1)),
    };
  }

  private updateColumnTags(column: DG.Column | undefined): void {
    if (column) {
      column.setTag(TAGS.LOW, column.getTag(TAGS.LOW) ?? DEFAULTS.LOW);
      column.setTag(TAGS.HIGH, column.getTag(TAGS.HIGH) ?? DEFAULTS.HIGH);
      column.setTag(TAGS.WEIGHT, column.getTag(TAGS.WEIGHT) ?? this.generateRandomNumber().toFixed(1));
    }
  }

  private initializeTreeGroups(): Map<string, DG.Column[]> {
    const groupMap = new Map<string, DG.Column[]>();

    this.columns.forEach(col => {
      const groupName = col.getTag(TAGS.GROUP_NAME);
      if (groupName) {
        const existingColumns = groupMap.get(groupName) ?? [];
        groupMap.set(groupName, [...existingColumns, col]);
      }
    });

    return groupMap;
  }

  public createTreeGroup(): HTMLElement {
    const inputs = ui.divV([]);
    const groupMap = this.initializeTreeGroups();
    const untaggedColumns = this.getUntaggedColumns();

    groupMap.forEach((groupColumns, groupName) => {
      const groupNode = this.tree.group(groupName);
      this.configureGroupNode(groupNode);

      groupColumns.forEach(col => {
        const colNode = groupNode.item(col.name);
        this.makeItemDraggable(colNode);
      });
    });

    this.tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
      if (node.parent.text !== '')
        this.updateInputs(inputs, node.text)
    });

    this.tree.onNodeContextMenu.subscribe((args: any) => {
      const menu: DG.Menu = args.args.menu;
      const node: DG.TreeViewNode = args.args.item;
      if (node instanceof DG.TreeViewGroup) {
        menu.item('Edit...', () => this.editGroup(node));
        menu.item('Delete', () => this.deleteGroup(node));
      }
      else
        menu.item('Add...', () => this.addGroup(node));
    });

    this.makeUntaggedColumnsDraggable(untaggedColumns);

    const generalInputs = ui.divV(generalProps.map(prop => this.createInputForProperty(prop, '')));
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
      subGroupProps.forEach(prop => {
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
        this.columns = this.columns.map(col => {
          if (col.getTag(TAGS.GROUP_NAME) === group.text)
            col.setTag(TAGS.GROUP_NAME, nameInput.stringValue);
          return col;
        });
  
        const sector = this.settings.sectors?.sectors.find(s => s.name === group.text);
        if (sector)
          sector.name = nameInput.stringValue;

        group.text = nameInput.stringValue;
      })
      .show();
  }

  private deleteGroup(group: DG.TreeViewGroup) {
    const groupChildren = group.children;
    group.remove();
    
    if (groupChildren) {
      groupChildren.forEach(child => {
      if (child instanceof DG.TreeViewNode) {
        const newChild = this.tree.item(child.text);
        delete this.columns.find((col) => col.name === child.text)?.tags[TAGS.GROUP_NAME];
        delete this.columns.find((col) => col.name === child.text)?.tags[TAGS.SECTOR_COLOR];
        this.makeItemDraggable(newChild);
      }
    });

    this.settings.sectors!.sectors = this.settings.sectors!.sectors.filter(sector => sector.name !== group.text);
    this.gc.grid.invalidate();
    }
  }

  private getUntaggedColumns(): DG.Column[] {
    return this.columns.filter(col => !col.getTag(TAGS.GROUP_NAME));
  }

  private makeUntaggedColumnsDraggable(columns: DG.Column[]): void {
    columns.forEach(col => {
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

export { VlaaiVisManager };