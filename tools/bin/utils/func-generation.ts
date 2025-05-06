/* eslint-disable no-unused-vars */
/* eslint-disable valid-jsdoc */
import { FuncMetadata, FuncParam } from './interfaces';

export const headerParams = ['name', 'description', 'tags', 'inputs', 'outputs'];

export enum pseudoParams {
  EXTENSION = 'extension',
  EXTENSIONS = 'extensions',
  INPUT_TYPE = 'inputType',
}

export enum FUNC_TYPES {
  APP = 'app',
  CELL_RENDERER = 'cellRenderer',
  FILE_EXPORTER = 'fileExporter',
  FILE_IMPORTER = 'file-handler',
  FILE_VIEWER = 'fileViewer',
  SETTINGS_EDITOR = 'packageSettingsEditor',
  VIEWER = 'viewer',
  FILTER = 'filter',
  AUTOSTART = 'autostart',
  INIT = 'init',
  EDITOR = 'editor',
  PANEL = 'panel',
  FOLDER_VIEWER = 'folderViewer',
  SEM_TYPE_DETECTOR = 'semTypeDetector',
  DASHBOARD = 'dashboard',
  FUNCTION_ANALYSIS = 'functionAnalysis',
  CONVERTER = 'converter'
}

export const typesToAnnotation : Record<string, string> = {
  'DataFrame': 'dataframe',
  'DG.DataFrame': 'dataframe',
  'Column': 'column',
  'DG.Column': 'column',
  'ColumnList': 'column_list',
  'DG.ColumnList': 'column_list',
  'FileInfo': 'file',
  'DG.FileInfo': 'file',
  'Uint8Array': 'blob',
  'number': 'double',
  'boolean': 'bool',
  'dayjs.Dayjs': 'datetime',
  'Dayjs': 'datetime',
  'graphics': 'graphics',
  'DG.View': 'view',
  'View': 'view',
  'DG.Widget': 'widget',
  'Widget': 'widget',
  'DG.FuncCall': 'funccall',
  'FuncCall': 'funccall',
  'DG.SemanticValue': 'semantic_value',
  'SemanticValue': 'semantic_value',
  'any': 'dynamic',
}

/** Generates an annotation header for a function based on provided metadata. */
export function getFuncAnnotation(data: FuncMetadata, comment: string = '//', sep: string = '\n'): string {
  const isFileViewer = data.tags?.includes(FUNC_TYPES.FILE_VIEWER) ?? false;
  const isFileImporter = data.tags?.includes(FUNC_TYPES.FILE_IMPORTER) ?? false;
  let s = '';
  if (data.name)
    s += `${comment}name: ${data.name}${sep}`;
  if (pseudoParams.EXTENSION in data && data.tags != null && data.tags.includes(FUNC_TYPES.FILE_EXPORTER))
    s += `${comment}description: Save as ${data[pseudoParams.EXTENSION]}${sep}`;
  else if (data.description)
    s += `${comment}description: ${data.description}${sep}`;
  if (data.tags) {
    s += `${comment}tags: ${isFileViewer && data[pseudoParams.EXTENSIONS] ?
      data.tags.concat(data[pseudoParams.EXTENSIONS].map((ext: string) => 'fileViewer-' + ext)).join(', ') :
      data.tags.join(', ')}${sep}`;
  }
  for(let input of data.inputs ?? [])
  {
    if(!input)
      continue;
    let type = input?.type;
    let isArray = false;
    if(type?.includes(`[]`)){
      type = type.replace(/\[\]$/, '');
      isArray = true;
    }
    const annotationType = typesToAnnotation[type ?? ''];
    if((input?.options as any)?.type)
      type = (input?.options as any)?.type;
    else if(annotationType)
    {
      if(isArray)
        type = `list<${annotationType}>`;
      else
        type = annotationType;
    }
    else
      type = 'dynamic';
    console.log(input);
    const options = ((input?.options as any)?.options? buildStringOfOptions((input.options as any).options ?? {}) : '');
    const functionName  = ((input.options as any)?.name ?  ((input?.options as any)?.name ?? ` ${input.name?.replaceAll('.', '')}`) : '')?.trim();
    s += comment + 'input: ' + type + ' ' + functionName + (input.defaultValue !== undefined ? `= ${input.defaultValue}` : '') + ' ' + options.replaceAll('"', '\'') + sep;
  }
  if (data.outputs) {
    for (const output of data.outputs)
      s += comment + 'output: ' + output.type + (output.name ? ` ${output.name}` : '') + sep;
  }
  
  if (data.meta) {
    for(let entry of Object.entries(data.meta))
      s += `${comment}meta.${entry[0]}: ${entry[1]}${sep}`;
  }

  for (const parameter in data) {
    if (parameter === pseudoParams.EXTENSION || parameter === pseudoParams.INPUT_TYPE || parameter === 'meta')
      continue;
    else if (parameter === pseudoParams.EXTENSIONS) {
      if (isFileViewer)
        continue;
      s += `${comment}meta.ext: ${data[parameter]}${sep}`;
    } else if (!headerParams.includes(parameter))
      s += `${comment}meta.${parameter}: ${data[parameter]}${sep}`;

  }
  return s;
}

function buildStringOfOptions(options: any){
  let optionsInString : string[] = [];
  for (const [key, value] of Object.entries(options ?? {})) {
    let val = value;
    if(Array.isArray(value))
      val = JSON.stringify(value);
    optionsInString.push(`${key}: ${val}`);
  }
  return `{ ${optionsInString.join('; ')} }`;
}

export const reservedDecorators: { [decorator: string]: { metadata: FuncMetadata, genFunc: Function } } = {
  viewer: {
    metadata: {
      tags: [FUNC_TYPES.VIEWER],
      inputs: [],
      outputs: [{ name: 'result', type: 'viewer' }],
    },
    genFunc: generateClassFunc,
  },
  filter: {
    metadata: {
      tags: [FUNC_TYPES.FILTER],
      inputs: [],
      outputs: [{ name: 'result', type: 'filter' }],
    },
    genFunc: generateClassFunc,
  },
  cellRenderer: {
    metadata: {
      tags: [FUNC_TYPES.CELL_RENDERER],
      inputs: [],
      outputs: [{ name: 'renderer', type: 'grid_cell_renderer' }],
    },
    genFunc: generateClassFunc,
  },
  fileExporter: {
    metadata: {
      tags: [FUNC_TYPES.FILE_EXPORTER],
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  fileImporter: {
    metadata: {
      tags: [FUNC_TYPES.FILE_IMPORTER],
      inputs: [{ name: 'content', type: 'string' }],
      outputs: [{ name: 'tables', type: 'list' }],
    },
    genFunc: generateFunc,
  },
  fileViewer: {
    metadata: {
      tags: [FUNC_TYPES.FILE_VIEWER],
      inputs: [{ name: 'f', type: 'file' }],
      outputs: [{ name: 'v', type: 'view' }],
    },
    genFunc: generateFunc,
  },
  settingsEditor: {
    metadata: {
      tags: [FUNC_TYPES.SETTINGS_EDITOR],
      inputs: [],
      outputs: [{ name: 'result', type: 'widget' }],
    },
    genFunc: generateFunc,
  },
  func: {
    metadata: {
      tags: [],
      inputs: [],
      outputs: [{ name: 'result', type: 'dynamic' }],
    },
    genFunc: generateFunc,
  },
  app: {
    metadata: {
      tags: [FUNC_TYPES.APP],
      inputs: [],
      outputs: [{ name: 'result', type: 'view' }],
    },
    genFunc: generateFunc,
  },
  autostart: {
    metadata: {
      tags: [FUNC_TYPES.AUTOSTART],
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  init: {
    metadata: {
      tags: [FUNC_TYPES.INIT],
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  editor: {
    metadata: {
      tags: [FUNC_TYPES.EDITOR],
      inputs: [{ name: 'call', type: 'funccall' }],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  panel: {
    metadata: {
      tags: [FUNC_TYPES.PANEL],
      inputs: [],
      outputs: [{ name: 'result', type: 'widget' }],
    },
    genFunc: generateFunc,
  },
  folderViewer: {
    metadata: {
      tags: [FUNC_TYPES.FOLDER_VIEWER],
      inputs: [{ name: 'folder', type: 'file' }, { name: 'files', type: 'list<file>' }],
      outputs: [{ name: 'result', type: 'widget' }],
    },
    genFunc: generateFunc,
  },
  semTypeDetector: {
    metadata: {
      tags: [FUNC_TYPES.SEM_TYPE_DETECTOR],
      inputs: [{ name: 'col', type: 'column' }],
      outputs: [{ name: 'result', type: 'string' }],
    },
    genFunc: generateFunc,
  },
  dashboard: {
    metadata: {
      tags: [FUNC_TYPES.DASHBOARD],
      inputs: [],
      outputs: [{ name: 'result', type: 'widget' }],
    },
    genFunc: generateFunc,
  },
  functionAnalysis: {
    metadata: {
      tags: [FUNC_TYPES.FUNCTION_ANALYSIS],
      inputs: [],
      outputs: [{ name: 'result', type: 'view' }],
    },
    genFunc: generateFunc,
  },
  converter: {
    metadata: {
      tags: [FUNC_TYPES.CONVERTER],
      inputs: [{ name: 'value', type: 'dynamic' }],
      outputs: [{ name: 'result', type: 'dynamic' }],
    },
    genFunc: generateFunc,
  },

};

/** Generates a DG function that instantiates a class. */
export function generateClassFunc(annotation: string, className: string, sep: string = '\n'): string {
  return annotation + `export function _${className}() {${sep}  return new ${className}();${sep}}${sep.repeat(2)}`;
}

/** Generates a DG function. */
export function generateFunc(annotation: string, funcName: string, sep: string = '\n', className: string = '', inputs: FuncParam[] = []): string 
{
  let funcSigNature = (inputs.map((e)=>`${e.name}: ${e.type}`)).join(', ');
  let funcArguments = (inputs.map((e)=>e.name)).join(', ');
  return annotation + `export function _${funcName}(${funcSigNature}) {${sep}  return ${className.length > 0 ? `${className}.` : ''}${funcName}(${funcArguments});${sep}}${sep.repeat(2)}`;
}

export function generateImport(className: string, path: string, sep: string = '\n'): string {
  return `import {${className}} from '${path}';${sep}`;
}

export function generateExport(className: string, sep: string = '\n'): string {
  return `export {${className}};${sep}`;
}
