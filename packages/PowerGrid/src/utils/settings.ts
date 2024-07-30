import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {PieChartSettings,Subsector} from '../sparklines/piechart';
import {generalProps, groupProps, subGroupProps} from './properties';

function generateRandomNumber(): number {
  return Math.random();
}

function getSectorProperty(settings: PieChartSettings, name: string, propertyName: string): any {
  for (const sector of settings.sectors!.sectors) {
    if (sector.name === name) {
      return sector[propertyName as keyof typeof sector];
    }
    for (const subsector of sector.subsectors) {
      if (subsector.name === name) {
        return subsector[propertyName as keyof typeof subsector];
      }
    }
  }
  return name === '' ? (settings.sectors as any)[propertyName] : undefined;
}

function updateSectorProperty(
  settings: PieChartSettings,
  name: string,
  propertyName: string,
  value: any,
  columns: DG.Column[] // Add tag parameter to match column tag with the name
): void {
  let sectorUpdated = false;
  let columnUpdated = false;

  // Create a map for fast column lookup by name
  const columnMap = new Map <string, DG.Column> ();
  columns.forEach(col => columnMap.set(col.name, col));

  // Update the property in the sector or subsector
  for (const sector of settings.sectors!.sectors) {
    if (sector.name === name) {
      sector[propertyName as keyof typeof sector] = value;
      sectorUpdated = true;
      break;
    }
    for (const subsector of sector.subsectors) {
      if (subsector.name === name) {
        (subsector[propertyName as keyof typeof subsector] as any) = value;
        sectorUpdated = true;
        break;
      }
    }
    if (sectorUpdated) break;
  }

  if (name === '') {
    (settings.sectors as any)[propertyName] = value;
  }

  // Update column tags using the map for direct access
  const column = columnMap.get(name);
  if (propertyName !== 'sectorColor') {
    if (column) {
      column.setTag(`.${propertyName}`, value);
      columnUpdated = true;
    } else {
      console.warn(`Column with name ${name} not found.`);
    }
  } else {
    const colsToUpdate = columns.map((c) => {
      if (c.getTag('.group-name') === name)
        return c;
    });
    for (const col of colsToUpdate)
      col?.setTag(`.${propertyName}`, value);
  }
}

function createInputForProperty(property: any, name: string, settings: PieChartSettings, gc: DG.GridColumn, columns: DG.Column[]): HTMLElement {
  const prop = DG.Property.fromOptions(property.property);
  const input = DG.InputBase.forProperty(prop, {});

  const value = getSectorProperty(settings, name, property.property.name);
  input.value = value !== undefined ? value : property.object[property.property.name];
  updateSectorProperty(settings, name, property.property.name, input.value, columns);

  input.onChanged(() => {
    updateSectorProperty(settings, name, property.property.name, input.value, columns);
    settings.sectors = settings.sectors; // Trigger change detection
    gc.grid.invalidate();
  });

  return input.root;
}

function makeItemDraggable(item: DG.TreeViewNode < any > ): void {
  ui.makeDraggable(item.root, {
    getDragObject: () => item,
    getDragCaption: () => `You are dragging ${item.text}`
  });
}

function makeGroupDroppable(groupNode: DG.TreeViewGroup, settings: PieChartSettings, columns: DG.Column[], tag: string, tree: DG.TreeViewGroup): void {
  ui.makeDroppable(groupNode.root, {
    acceptDrop: () => true,
    doDrop: (draggedItem: any, _) => {
      const itemText = draggedItem.text;
      // Find and remove the item from its original group
      // Remove the item from its original group
      tree.items.forEach((group) => {
        if (group.text === itemText) {
          group.remove();
          const originalGroup = settings.sectors?.sectors.find(sector => sector.name === group.parent.text);
          if (originalGroup)
            originalGroup.subsectors = originalGroup.subsectors.filter(subsector => subsector.name !== itemText);
        }
      });

      // Add the item to the new group
      let newGroup = settings.sectors?.sectors.find(sector => sector.name === groupNode.text);
      const column = columns.find(col => col.name === itemText);

      if (!newGroup) {
        settings.sectors?.sectors.push({
          name: groupNode.text,
          sectorColor: column?.getTag('.sectorColor') || defaultGroupProps["sectorColor"],
          subsectors: []
        });
      }

      newGroup = settings.sectors?.sectors.find(sector => sector.name === groupNode.text);
      column?.setTag('.low', column?.getTag('.low') || '0');
      column?.setTag('.high', column?.getTag('.high') || '1');
      column?.setTag('.weight', column?.getTag('.weight') || generateRandomNumber().toFixed(1));
      if (!newGroup?.subsectors.some(subsector => subsector.name === itemText)) {
        newGroup?.subsectors.push({
          name: itemText,
          lowThreshold: parseFloat(column?.getTag('.low') !),
          highThreshold: parseFloat(column?.getTag('.high') !),
          weight: parseFloat(column?.getTag('.weight') !),
        });
      }

      if (column) {
        column.setTag(tag, groupNode.text);
      }

      // Make the new item draggable
      const newItem = groupNode.item(itemText);
      makeItemDraggable(newItem);
    }
  });
}

function initializeTreeGroups(settings: PieChartSettings, gc: DG.GridColumn, columns: DG.Column[], tag: string): Map < string, DG.Column[] > {
  const groupMap = new Map < string,
    DG.Column[] > ();
  const propertiesMap = new Map < string,
    any > ();

  columns.forEach(col => {
    const groupName = col.getTag(tag);
    if (groupName) {
      if (!groupMap.has(groupName)) {
        groupMap.set(groupName, []);
        propertiesMap.set(groupName, {});
      }
      groupMap.get(groupName) !.push(col);
      propertiesMap.set(col.name, {});
    }
  });

  return groupMap;
}

export function createTreeGroup(settings: PieChartSettings, gc: DG.GridColumn): HTMLElement {
  const df = grok.shell.tv.dataFrame;
  const tag = '.group-name';
  const tree = ui.tree();
  const inputs = ui.divV([]);
  const columns = df.columns.byNames(settings.columnNames);

  // Initialize tree groups and map
  const groupMap = initializeTreeGroups(settings, gc, columns, tag);
  const untaggedColumns: string[] = [];

  // Identify columns without tags
  columns.forEach(col => {
    if (!col.getTag(tag)) {
      untaggedColumns.push(col.name);
    }
  });

  // Set up tree groups
  groupMap.forEach((groupColumns, groupName) => {
    const groupNode = tree.group(groupName);
    const colorPicker = createInputForProperty(groupProps[0], groupName, settings, gc, columns).getElementsByClassName('ui-input-options')[0];
    //(colorPicker as HTMLElement).style.margin = '5px!important';
    (colorPicker as HTMLElement).style.cssText += ('margin: 5px!important');
    //const colorPicker = ui.input.color('').root;
    const referenceNode = groupNode.root.getElementsByClassName('d4-tree-view-group-label')[0];
    referenceNode.insertAdjacentElement('beforebegin', colorPicker);

    groupNode.expanded = false;
    makeGroupDroppable(groupNode, settings, columns, tag, tree);

    groupColumns.forEach(col => {
      const colNode = groupNode.item(col.name);
      colNode.checked = false;
      makeItemDraggable(colNode);
    });
  });

  // Handle node selection changes
  tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
    ui.empty(inputs);
    const isSubGroup = settings.columnNames.includes(node.text);
    if (!isSubGroup)
      return;
    subGroupProps.forEach(prop => {
      inputs.appendChild(createInputForProperty(prop, node.text, settings, gc, columns));
    });
  });

  // Create general inputs
  const generalInputs = ui.divV(generalProps.map(prop => createInputForProperty(prop, '', settings, gc, columns)));

  // Make untagged columns draggable
  untaggedColumns.forEach(colName => {
    const colNode = tree.item(colName);
    makeItemDraggable(colNode);
  });

  tree.root.style.overflow = 'hidden';
  //@ts-ignore
  tree.currentItem = tree.children[0];
  generalInputs.style.marginRight = '10px';
  const resultingDiv = ui.divV([generalInputs, ui.divH([tree.root, inputs], 'ui-form')]);
  resultingDiv.style.marginTop = '10px';

  const addGroup = () => {
    const nameInput = ui.input.string('Name');
    ui.dialog('Add group')
      .add(nameInput)
      .onOK(() => {
        const lastGroupIndex = findLastGroupIndex(tree.children);
        const newGroup = tree.group(nameInput.value, null, true, lastGroupIndex);
        const colorPicker = createInputForProperty(groupProps[0], nameInput.value, settings, gc, columns).getElementsByClassName('ui-input-options')[0];
        (colorPicker as HTMLElement).style.cssText += ('margin: 5px!important');
        const referenceNode = newGroup.root.getElementsByClassName('d4-tree-view-group-label')[0];
        referenceNode.insertAdjacentElement('beforebegin', colorPicker);
        newGroup.expanded = false;
        makeGroupDroppable(newGroup, settings, columns, tag, tree);
      })
      .show();
  }

  const deleteGroup = () => {
    const nameInput = ui.input.string('Name');

    ui.dialog('Delete group')
      .add(nameInput)
      .onOK(() => {
        const groupName = nameInput.value;
        const groupToDelete = tree.children.find(node =>
          node instanceof DG.TreeViewGroup && node.text === groupName
        ) as DG.TreeViewGroup;

        const groupChildren = groupToDelete.children;
        groupToDelete.remove();

        if (groupChildren) {
          groupChildren.forEach(child => {
            if (child instanceof DG.TreeViewNode) {
              //tree.root.appendChild(child.root);
              //tree.items.push(child);
              const newChild = tree.item(child.text);
              console.log('tree items');
              console.log(tree.items);
              delete columns.find((col) => col.name === child.text)?.tags['.group-name'];
              delete columns.find((col) => col.name === child.text)?.tags['.sectorColor'];
              //tree.items.push(child);
              makeItemDraggable(newChild);
            }
          });

          settings.sectors!.sectors = settings.sectors!.sectors.filter(sector => sector.name !== groupName);
          gc.grid.invalidate();
        }
      })
      .show();
  }

  resultingDiv.oncontextmenu = (e) => {
    DG.Menu.popup()
      .item('Add group', () => addGroup())
      .item('Delete group', () => deleteGroup())
      .show();
    e.preventDefault();
  };
  return resultingDiv;
}

function findLastGroupIndex(children: DG.TreeViewNode < any > []): number {
  for (let i = children.length - 1; i >= 0; i--) {
    if (children[i] instanceof DG.TreeViewGroup) {
      return i + 1;
    }
  }
  return 0; // No group found
}

const defaultGeneralProps = generalProps.reduce((acc, prop) => {
    acc[prop.property.name] = (prop.object as any)[prop.property.name];
    return acc;
  }, {} as Record <string, any>);

const defaultGroupProps = groupProps.reduce((acc, prop) => {
    acc[prop.property.name] = (prop.object as any)[prop.property.name];
    return acc;
  }, {} as Record <string, any>);

export function initializeSectors(columnNames: string[], dataFrame: DG.DataFrame): {
  lowerBound: number;
  upperBound: number;
  sectors: {
    name: string;
    sectorColor: string;
    subsectors: Subsector[];
  } [];
  values: string | null;
} {
  const groupMap = new Map < string,
    DG.Column[] > ();
  const columns = dataFrame.columns.byNames(columnNames);
  const tag = '.group-name';

  columns.forEach(col => {
    const groupName = col.getTag(tag);
    if (groupName) {
      if (!groupMap.has(groupName)) {
        groupMap.set(groupName, []);
      }
      groupMap.get(groupName) !.push(col);

      // Set tags for the column if not already set
      col.setTag('.low', col.getTag('.low') || '0');
      col.setTag('.high', col.getTag('.high') || '1');
      col.setTag('.weight', col.getTag('.weight') || generateRandomNumber().toFixed(1));
    }
  });

  const sectors = Array.from(groupMap.entries()).map(([groupName, columns]) => ({
    name: groupName,
    sectorColor: columns[0].getTag('.sectorColor') || defaultGroupProps["sectorColor"],
    subsectors: columns.map(col => ({
      name: col.name,
      lowThreshold: parseFloat(col.getTag('.low')) || 0,
      highThreshold: parseFloat(col.getTag('.high')) || 1,
      weight: parseFloat(col.getTag('.weight')) || parseFloat(generateRandomNumber().toFixed(1)),
    }))
  }));

  return {
    lowerBound: defaultGeneralProps["lowerBound"] || DEFAULT_LOWER_VALUE,
    upperBound: defaultGeneralProps["upperBound"] || DEFAULT_UPPER_VALUE,
    sectors,
    values: null
  };
}

const DEFAULT_LOWER_VALUE = 0.7;
const DEFAULT_UPPER_VALUE = 0.9;