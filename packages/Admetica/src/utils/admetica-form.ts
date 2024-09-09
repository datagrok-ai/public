import * as ui from 'datagrok-api/ui';
import { tablePieChartIndexMap } from './admetica-utils';

// Function to create an element state for each column
function createElementState(table: string, left: number, top: number, column: string) {
  return [
    {
      "left": left + getTextWidth(column),  // Position the field to the right of the HTML input
      "top": top,
      "width": 100,
      "height": 20,
      "type": "field",
      "viewerSettings": {
        "table": table,
        "column": column,
        "format": null
      }
    },
    {
      "left": left,
      "top": top,
      "width": getTextWidth(column),
      "height": 20,
      "type": "html",
      "viewerSettings": {
        "markup": "<input type=\"text\" class=\"d4-sketch-column-name ui-input-editor\">",
        "table": table,
        "column": column,
        "input-value": column
      }
    }
  ];
}

function getTextWidth(text: string) {
  const canvas = ui.canvas();
  const context = canvas.getContext('2d');
  context!.font = '13px Roboto, sans-serif';
  const width = context!.measureText(text).width;

  return width + 26;
}
  
// Function to create a header element for a category
function createCategoryHeader(left: number, top: number, categoryName: string) {
  return {
    "left": left,
    "top": top,
    "width": 280,
    "height": 20, // Increased height to ensure adequate spacing
    "type": "html",
    "viewerSettings": {
      "markup": "<input type=\"text\" class=\"d4-sketch-column-name ui-input-editor\">",
      "input-value": `                         ${categoryName}`,
      "backgroundColor": 4293717745
    }
  };
}
  
// Function to generate the form state with dynamic vertical placement
export function generateFormState(table: string, categories: { [s: string]: string[]; }, molColName: string, addPiechart: boolean, maxHeight: number = 600) {
  const elementStates = [];
  const leftOffset = 10; // Starting left offset
  const colWidth = 300; // Width allocated for each category set
  const headerHeight = 30; // Space for headers
  const rowHeight = 30; // Height of each row including spacing
  const columnGap = 70; // Larger gap between columns
  
  let currentTopOffset = 180;  // Starting top offset, under "smiles" and "sparkline-cell"
  let currentLeftOffset = leftOffset;
  let piechartIndex = tablePieChartIndexMap.get(table);

  // Add the predefined "smiles" and "sparkline-cell" elements at their specific positions
  elementStates.push(
    {
      "left": 3,
      "top": 21,
      "width": 183,
      "height": 139,
      "type": "field",
      "viewerSettings": {
        "table": table,
        "column": molColName,
        "format": null
      }
    }
  );

  if (addPiechart)
    elementStates.push({
      "left": 154,
      "top": 22,
      "width": 209,
      "height": 140,
      "type": "sparkline-cell",
      "viewerSettings": {
        "table": table,
        "column": piechartIndex === 1 ? "piechart" : `piechart (${piechartIndex})`
      }
    });
  
  // Generate other elements below the "smiles" and "sparkline-cell"
  for (const [category, columns] of Object.entries(categories)) {
    // Calculate the height of the category (header + rows)
    const categoryHeight = headerHeight + columns.length * rowHeight;
  
    // If adding this category exceeds the max height, start a new column
    if (currentTopOffset + categoryHeight > maxHeight) {
      currentLeftOffset += colWidth + columnGap;  // Move to the next column with an additional gap
      currentTopOffset = 30;  // Reset to the top offset for the new column
    }
  
    // Add a header for the category
    const header = createCategoryHeader(currentLeftOffset, currentTopOffset, category);
    elementStates.push(header);
    currentTopOffset += headerHeight; // Move down to place elements under the header
  
    // Add elements for the current category's columns
    columns.forEach((column: string, index: number) => {
      const elementState = createElementState(table, currentLeftOffset, currentTopOffset + index * rowHeight, column);
      elementStates.push(...elementState); // Add both HTML and field elements
    });
  
    currentTopOffset += columns.length * rowHeight + 10; // Update top offset for the next category, adding a little extra spacing
  }
  
  return {
    "#type": "SketchState",
    "elementStates": elementStates
  };
}
  
export function createCategoryModelMapping(properties: any, updatedModelNames: string[]) {
  const categoryModelMapping = {};
  
  // Iterate over each subgroup in the properties object
  properties.subgroup.forEach((subgroup: { name: string; models: { name: string }[]; }) => {
    const categoryName = subgroup.name;  // Get the category name from the subgroup
  
    // Filter the updatedModelNames to only those that are present in subgroup.models
    const modelNames = updatedModelNames
      .filter(name => subgroup.models.some(model => name.includes(model.name)));
  
    // Add the category and its corresponding model names to the mapping
    if (modelNames.length > 0) {
      (categoryModelMapping as any)[categoryName] = modelNames;
    }
  });  
  
  return categoryModelMapping;
}