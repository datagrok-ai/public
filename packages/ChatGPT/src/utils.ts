/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

const AI_SQL_QUERY_ABORT_EVENT = 'd4-ai-generation-abort';
const AI_PANEL_TOGGLE_EVENT = 'd4-ai-panel-toggle';

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

  viewer.getProperties().forEach((prop) => {
    properties[prop.name] = {
      type: prop.type,
      description: prop.description
    };
  });

  return {
    viewerType,
    properties
  };
}

export function getDataFrameDescription(dataFrame: DG.DataFrame): IDataFrameDescription {
  const columns: { [columnName: string]: IProp } = {};

  for (const column of dataFrame.columns.all) {
    columns[column.name] = {type: column.type};

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

    viewer.getProperties().forEach((prop) => {
      properties[prop.name] = {
        type: prop.type,
        description: prop.description
      };
    });

    descriptions[viewerType] = properties;
  }

  return descriptions;
}

const commonProps = ['rowSource', 'filter', 'table'];
/**
 * Returns a list of viewer types, along with the detailed descriptions of their properties.
 * This is used to generate the AI assistant's context.
 */
export function getViewerDescriptionsString(): string {
  if (viewerDescriptions !== '')
    return viewerDescriptions;


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
    } catch (e) {
      grok.shell.error('Error getting viewer descriptions: ' + e);
    }

    result += '\n';
  }

  viewerDescriptions = result;
  return result;
}

export function getCurrentViewersString(view: DG.TableView): string {
  let result = `Current viewers in the view:\n`;
  for (const viewer of view.viewers) {
    result += `${viewer.type} with properties: \n`;
    for (const prop of viewer.getProperties()) {
      if (!commonProps.includes(prop.name) && (prop.category == 'Data' ||
        prop.name.endsWith('ColumnName')) && viewer.props[prop.name]
      )
        result += `${prop.name}: ${viewer.props[prop.name]?.toString()}` + '\n';
    }
    result += '\n\n';
  }
  return result;
}

export type AbortPointer = {
  aborted: boolean;
}

export function dartLike<T extends any>(obj: T) {
  return {
    set: function<K extends keyof T>(key: K, value: T[K]) {
      (obj as any)[key] = value;
      return this;
    },
    value: obj as T,
  };
}

export function getAIAbortSubscription() {
  return rxjs.merge(grok.events.onEvent(AI_SQL_QUERY_ABORT_EVENT), grok.events.onCustomEvent(AI_SQL_QUERY_ABORT_EVENT));
}

export function getAIPanelToggleSubscription() {
  return rxjs.merge(grok.events.onEvent(AI_PANEL_TOGGLE_EVENT), grok.events.onCustomEvent(AI_PANEL_TOGGLE_EVENT));
}

export function fireAIAbortEvent() {
  grok.events.fireCustomEvent(AI_SQL_QUERY_ABORT_EVENT, null);
}

export function fireAIPanelToggleEvent(v: DG.View | DG.ViewBase) {
  grok.events.fireCustomEvent(AI_PANEL_TOGGLE_EVENT, v);
}

export function findLast<T, K extends T>(array: T[], predicate: (value: T, index: number, obj: T[]) => value is K): K | undefined {
  for (let i = array.length - 1; i >= 0; i--) {
    const value = array[i];
    if (predicate(value, i, array))
      return value;
  }
  return undefined;
}
