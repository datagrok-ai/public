import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import { PieChartSettings, Sector, Subsector } from '../sparklines/piechart';
import { generalProps, groupProps, subGroupProps } from './properties';

enum SectorType {
  SECTOR = 'sector',
  SUBSECTOR = 'subsector'
}

const Constants = {
  TAG_PREFIX: '.',
  GROUP_NAME_TAG: '.group-name', //move to tags
  SECTOR_COLOR_PROPERTY: 'sectorColor'
};

const TAGS = {
  SECTOR_COLOR: '.sectorColor',
  LOW: '.low',
  HIGH: '.high',
  WEIGHT: '.weight',
  GROUP_NAME: '.group-name'
};

const DEFAULTS = {
  LOW: '0',
  HIGH: '1',
  WEIGHT: '0.0'
};

const defaultGeneralProps = generalProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record <string, any>);

const defaultGroupProps = groupProps.reduce((acc, prop) => {
  acc[prop.property.name] = (prop.object as any)[prop.property.name];
  return acc;
}, {} as Record <string, any>);

class PieChartManager {
  private settings: PieChartSettings;
  private gc: DG.GridColumn;
  private columns: DG.Column[];
  private tree: DG.TreeViewGroup;
  private readonly tag: string = '.group-name';

  constructor(settings: PieChartSettings, gc: DG.GridColumn) {
    this.settings = settings;
    this.gc = gc;
    this.columns = grok.shell.tv.dataFrame.columns.byNames(settings.columnNames);
    this.tree = ui.tree();

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
      const groupName = col.getTag(this.tag);
      if (groupName) {
        const existingColumns = groupMap.get(groupName) ?? [];
        groupMap.set(groupName, [...existingColumns, col]);
        this.updateColumnTags(col);
      }
    });

    const sectors = Array.from(groupMap.entries()).map(([groupName, columns]) => ({
      name: groupName,
      sectorColor: columns[0].getTag(TAGS.SECTOR_COLOR) || defaultGroupProps["sectorColor"],
      subsectors: columns.map(col => ({
        name: col.name,
        lowThreshold: parseFloat(col.getTag(TAGS.LOW)) || 0,
        highThreshold: parseFloat(col.getTag(TAGS.HIGH)) || 1,
        weight: parseFloat(col.getTag(TAGS.WEIGHT)) || parseFloat(this.generateRandomNumber().toFixed(1)),
      }))
    }));

    return {
      lowerBound: defaultGeneralProps["lowerBound"],
      upperBound: defaultGeneralProps["upperBound"],
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
      propertyName === Constants.SECTOR_COLOR_PROPERTY 
        ? this.updateSectorColorTags(name, value)
        : this.updateColumnTag(name, propertyName, value, columnMap);
    } else if (name === '')
      (this.settings.sectors as any)[propertyName] = value;
  }

  private updateSectorColorTags(name: string, value: any): void {
    this.columns
      .filter(col => col.getTag(Constants.GROUP_NAME_TAG) === name)
      .forEach(col => col.setTag(`${Constants.TAG_PREFIX}${Constants.SECTOR_COLOR_PROPERTY}`, value));
  }

  private updateColumnTag(name: string, propertyName: keyof (Sector | Subsector), value: any, columnMap: Map<string, DG.Column>): void {
    const column = columnMap.get(name);
    column?.setTag(`${Constants.TAG_PREFIX}${propertyName}`, value);
  }

  private createInputForProperty(property: any, name: string): HTMLElement {
    const prop = DG.Property.fromOptions(property.property);
    const input = DG.InputBase.forProperty(prop, {});
    const value = this.getSectorProperty(name, property.property.name as keyof (Sector | Subsector));
    input.value = value !== undefined ? value : property.object[property.property.name];
    
    input.onChanged(() => {
      this.updateSectorProperty(name, property.property.name as keyof (Sector | Subsector), input.value);
      this.gc.grid.invalidate();
    });

    return input.root;
  }

  private makeItemDraggable(item: DG.TreeViewNode<any>): void {
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

    if (!newGroup.subsectors.some(subsector => subsector.name === itemText))
      newGroup.subsectors = [...newGroup.subsectors, this.createSubsector(column)];

    if (column)
      column.setTag(TAGS.GROUP_NAME, groupNode.text);

    const newItem = groupNode.item(itemText);
    this.makeItemDraggable(newItem);
  }

  private getOrCreateNewGroup(groupName: string): Sector {
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
      lowThreshold: parseFloat(column?.getTag(TAGS.LOW) ?? DEFAULTS.LOW),
      highThreshold: parseFloat(column?.getTag(TAGS.HIGH) ?? DEFAULTS.HIGH),
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
      const groupName = col.getTag(this.tag);
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

    this.tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => this.updateInputs(inputs, node.text));

    const generalInputs = this.createGeneralInputs();
    this.makeUntaggedColumnsDraggable(untaggedColumns);

    const resultingDiv = this.createResultingDiv(inputs, generalInputs);
    resultingDiv.appendChild(this.createAddGroupButton());

    return resultingDiv;
  }

  private configureGroupNode(groupNode: DG.TreeViewGroup): void {
    const colorPicker = this.createColorPicker(groupNode.text);
    const referenceNode = groupNode.root.getElementsByClassName('d4-tree-view-group-label')[0];
    referenceNode.insertAdjacentElement('beforebegin', colorPicker);

    groupNode.expanded = false;
    this.makeGroupDroppable(groupNode);
  }

  private createColorPicker(groupName: string): HTMLElement {
    const colorPicker = this.createInputForProperty(groupProps[0], groupName)
      .getElementsByClassName('ui-input-options')[0] as HTMLElement;
    colorPicker.style.margin = '5px!important';
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

  private createGeneralInputs(): HTMLElement {
    const container = ui.divV([this.tree.root, generalProps.map(prop => this.createInputForProperty(prop, ''))]);
    container.style.marginRight = '10px';
    return container;
  }

  private createResultingDiv(inputs: HTMLElement, generalInputs: HTMLElement): HTMLElement {
    const container = ui.divH([generalInputs, inputs]);
    container.style.marginTop = '10px';
    return container;
  }

  private createAddGroupButton(): HTMLElement {
    const addGroup = () => {
      const nameInput = ui.input.string('Name');
      ui.dialog('Add group')
        .add(nameInput)
        .onOK(() => {
          const newGroup = this.tree.group(nameInput.value, null, false, this.findLastGroupIndex());
          this.configureGroupNode(newGroup);
        })
        .show();
    };

    return ui.button('Add Group', addGroup);
  }

  private getUntaggedColumns(): DG.Column[] {
    return this.columns.filter(col => !col.getTag(this.tag));
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
    return this.tree.items.reduce((index, item, i) => item instanceof DG.TreeViewGroup ? Math.max(index, i) : index, 0);
  }
}

export { PieChartManager };