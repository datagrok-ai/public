import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const dummy = grok.data.testData('demog');
let viewerDescriptions: string = '';

type IProp = {
  type: string;
  description?: string;
}

type IViewerDescription = {
  viewerType: string;
  properties: { [propName: string]: IProp };
}

type IDataFrameDescription = {
  name: string;
  columns: { [columnName: string]: IProp };
}

export function getViewerDescription(viewerType: DG.ViewerType): IViewerDescription {
  const viewer = DG.Viewer.fromType(viewerType, dummy);
  const properties: { [propName: string]: IProp } = {};
  
  viewer.getProperties().forEach(prop => {
    properties[prop.name] = {
      type: prop.type,
      description: prop.description
    };
  });
  
  return {
    viewerType,
    properties
  }
}

export function getDataFrameDescription(dataFrame: DG.DataFrame): IDataFrameDescription {
  const columns: { [columnName: string]: IProp } = {};

  for (const column of dataFrame.columns.all) {
    columns[column.name] = { type: column.type };

    if (column.meta.description)
      columns[column.name].description = column.meta.description;
  }

  return {
    name: dataFrame.name,
    columns
  };
}

/** 
 * Returns a list of viewer types, along with the detailed descriptions of their properties.
 * This is used to generate the AI assistant's context.
 */
export function getViewerDescriptions(): {[viewerType: string]: { [propName: string]: IProp }} {
  grok.shell.info('gvd');
  const descriptions: {[viewerType: string]: { [propName: string]: IProp }} = {};
  
  for (const viewerType of DG.Viewer.CORE_VIEWER_TYPES) {
    const viewer = DG.Viewer.fromType(viewerType, dummy);
    const properties: { [propName: string]: IProp } = {};
    
    viewer.getProperties().forEach(prop => {
      properties[prop.name] = {
        type: prop.type,
        description: prop.description
      };
    });
    
    descriptions[viewerType] = properties;
  }
  
  return descriptions;
}

/** 
 * Returns a list of viewer types, along with the detailed descriptions of their properties.
 * This is used to generate the AI assistant's context.
 */
export function getViewerDescriptionsString(): string {
  if (viewerDescriptions !== '') 
    return viewerDescriptions;

  const commonProps = ['rowSource', 'filter', 'table'];
  let result = `Properties common for all viewers:
rowSource: string  // choices: "All", "Selected", "Filtered"
filter: string  // Formulas like "\${AGE} < 40" 

`;

  for (const viewerType of DG.Viewer.CORE_VIEWER_TYPES) {

    try {
      const viewer = DG.Viewer.fromType(viewerType, dummy);

      result += `${viewerType}\n`;
      for (const prop of viewer.getProperties()) {
        if (!commonProps.includes(prop.name) && (prop.category == 'Data' || prop.name.endsWith('ColumnName')))
          result += `${prop.name} ${prop.type}` + '\n';
      }
    }
    catch (e) {
      grok.shell.error('Error getting viewer descriptions: ' + e);
    }

    result += '\n';
  }

  viewerDescriptions = result;
  return result;
}