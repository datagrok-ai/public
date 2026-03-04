/* eslint-disable no-unused-vars */
/* eslint-disable valid-jsdoc */
import {FuncMetadata, FuncParam} from './interfaces';

export const headerParams = ['name', 'description', 'tags', 'inputs', 'outputs'];

export enum pseudoParams {
  EXTENSION = 'extension',
  EXTENSIONS = 'extensions',
  INPUT_TYPE = 'inputType',
}

const nonMetaData = [
  'sidebar',
  'editor',
  'friendlyName',
  'helpUrl', 
  'help-url', 
  'condition',
  'connection',
  'top-menu',
  'cache',
  'cache.invalidateOn',
  'test',
];

export const decoratorOptionToAnnotation = new Map<string, string>([
  ['helpUrl', 'help-url'],
  ['topMenu', 'top-menu'],
  ['metaUrl', 'meta.url'],
  ['optionsType', 'type'],
  ['cacheInvalidateOn', 'cache.invalidateOn'],
  ['scriptHandlerLanguage', 'scriptHandler.language'],
  ['scriptHandlerExtensions', 'scriptHandler.extensions'],
  ['scriptHandlerCommentStart', 'scriptHandler.commentStart'],
  ['scriptHandlerTemplateScript', 'scriptHandler.templateScript'],
  ['scriptHandlerCodeEditorMode', 'scriptHandler.codeEditorMode'],
  ['scriptHandlerVectorizationFunction', 'scriptHandler.vectorizationFunction'],
]);

export const dgAnnotationTypes: Record<string, string> = {
  INT: 'int',
  BIG_INT: 'bigint',
  FLOAT: 'double',
  NUM: 'num',
  QNUM: 'qnum',
  BOOL: 'bool',
  STRING: 'string',
  STRING_LIST: 'string_list',
  DATE_TIME: 'datetime',
  OBJECT: 'object',
  BYTE_ARRAY: 'byte_array',
  DATA_FRAME: 'dataframe',
  DATA_FRAME_LIST: 'dataframe_list',
  CELL: 'cell',
  COLUMN: 'column',
  COLUMN_LIST: 'column_list',
  GRAPHICS: 'graphics',
  FILE: 'file',
  BLOB: 'blob',
  ROW_FILTER: 'tablerowfiltercall',
  COLUMN_FILTER: 'colfiltercall',
  BIT_SET: 'bitset',
  MAP: 'map',
  DYNAMIC: 'dynamic',
  VIEWER: 'viewer',
  LIST: 'list',
  SEM_VALUE: 'semantic_value',
  FUNC: 'func',
  FUNC_CALL: 'funccall',
  PROPERTY: 'property',
  CATEGORICAL: 'categorical',
  NUMERICAL: 'numerical',
  GRID_CELL_RENDER_ARGS: 'GridCellRenderArgs',
  ELEMENT: 'element',
  VIEW: 'view',
  TABLE_VIEW: 'TableView',
  USER: 'User',
  MENU: 'Menu',
  PROJECT: 'Project',
  SEMANTIC_VALUE: 'semantic_value',
  EVENT_DATA: 'event_data',
  PROGRESS_INDICATOR: 'progressindicator',
  CREDENTIALS: 'Credentials',
  SCRIPT_ENVIRONMENT: 'ScriptEnvironment',
  NOTEBOOK: 'Notebook',
};

export enum FUNC_TYPES {
  APP = 'app',
  CELL_RENDERER = 'cellRenderer',
  FILE_EXPORTER = 'fileExporter',
  FILE_HANDLER = 'file-handler',
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
  CONVERTER = 'converter',
  MODEL = 'model',
  APP_TREE_BROWSER = 'appTreeBrowser'
}

export const typesToAnnotation: Record<string, string> = {
  'DataFrame': 'dataframe',
  'DG.DataFrame': 'dataframe',
  'Column': 'column',
  'DG.Column': 'column',
  'Column<any>': 'column',
  'DG.Column<any>': 'column',
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
  'DG.ViewBase': 'view',
  'View': 'view',
  'ViewBase': 'view',
  'DG.Widget': 'widget',
  'Widget': 'widget',
  'DG.FuncCall': 'funccall',
  'FuncCall': 'funccall',
  'DG.SemanticValue': 'semantic_value',
  'SemanticValue': 'semantic_value',
  'any': 'dynamic',
  'void': 'void',
  'string': 'string',
};

export const typesToAny: string[] = [
  'dayjs.Dayjs',
  'Dayjs',
];

/** Generates an annotation header for a function based on provided metadata. */
export function getFuncAnnotation(data: FuncMetadata, comment: string = '//', sep: string = '\n'): string {
  const isFileViewer = data.tags?.includes(FUNC_TYPES.FILE_VIEWER) ?? false;
  const isFileImporter = data.tags?.includes(FUNC_TYPES.FILE_HANDLER) ?? false;
  let s = '';
  if (data.name)
    s += `${comment}name: ${data.name}${sep}`;
  if (pseudoParams.EXTENSION in data && data.tags != null && data.tags.includes(FUNC_TYPES.FILE_EXPORTER))
    s += `${comment}description: Save as ${data[pseudoParams.EXTENSION]}${sep}`;
  else if (data.description)
    s += `${comment}description: ${data.description}${sep}`;
  if (data.tags && data.tags?.length > 0) {
    s += `${comment}tags: ${isFileViewer && data[pseudoParams.EXTENSIONS] ?
      data.tags.concat(data[pseudoParams.EXTENSIONS].map((ext: string) => 'fileViewer-' + ext)).join(', ') :
      data.tags.join(', ')}${sep}`;
  }

  for (const input of data.inputs ?? []) {
    if (!input)
      continue;
    let type = input?.type?.replaceAll(' ', '');
    if (type?.includes(`|undefined`)) {
      type = type.replaceAll(`|undefined`, '');
      // modify original payload
      input.optional = true;
      input.type = type;
    }
    let isArray = false;
    if (type?.includes(`[]`)) {
      type = type.replace(/\[\]$/, '');
      isArray = true;
    }

    const annotationType = typesToAnnotation[type ?? ''];
    if ((input?.options as any)?.type || (input?.options as any)?.options?.type)
      type = (input?.options as any)?.type ?? (input?.options as any)?.options?.type;
    else if (annotationType) {
      if (isArray)
        type = `list<${annotationType}>`;
      else
        type = annotationType;
    } else
      type = 'dynamic'; 
    
    // if ((input as any)?.options?.type === 'categorical' || (input as any)?.options?.options?.type === 'categorical')
    //   console.log(input);
    // console.log(input);

    const options = ((input?.options as any)?.options ? buildStringOfOptions((input as any).options ?? {}) : '');
    const functionName = ((input.options as any)?.name ? (input?.options as any)?.name : ` ${input.name?.replaceAll('.', '')}`)?.trim();
    
    // eslint-disable-next-line max-len
    s += comment + 'input: ' + type + ' ' + functionName + (input.defaultValue !== undefined ? `= ${input.defaultValue}` : '') + ' ' + options + sep;
  }
  if (data.outputs) {
    for (const output of data.outputs) {
      if (output.type !== 'void') {
      // eslint-disable-next-line max-len
        s += comment + 'output: ' + output.type + (output.name ? ` ${output.name}${output.options ? ` ${buildStringOfOptions(output)}` : ''}` : '') + sep;
      }
    }
  }

  if (data.meta) {
    for (const entry of Object.entries(data.meta))
      s += `${comment}meta.${decoratorOptionToAnnotation.get(entry[0]) ?? entry[0]}: ${entry[1]}${sep}`;
  }

  for (const parameter in data) {
    // eslint-disable-next-line max-len
    if (parameter === pseudoParams.EXTENSION || parameter === pseudoParams.INPUT_TYPE || parameter === 'meta' || parameter === 'isAsync' || parameter === 'test' || parameter === 'actualTypeObj')
      continue;
    else if (parameter === pseudoParams.EXTENSIONS) {
      if (isFileViewer)
        continue;
      s += `${comment}meta.ext: ${data[parameter]}${sep}`;
    } else if (!headerParams.includes(parameter)) {
      if (nonMetaData.includes(parameter))
        s += `${comment}${parameter}: ${data[parameter]}${sep}`;
      else
        s += `${comment}meta.${parameter}: ${data[parameter]}${sep}`;
    }
  }

  if (data.test) {
    if (typeof(data.test) === 'string') 
      s += `${comment}test: ${data.test}${sep}`;
    else {
      for (const entry of Object.entries(data.test)) {
        if (entry[0] === 'test' || entry[0] === 'wait')
          s += `${comment}`;
        else
          s += `, `;
        s += `${entry[0]}: ${entry[1]} `;
      }
      s += `${sep}`;
    }
  }
  return s;
}

export const inputOptionsNames = [
  'semType',
  'category',
  'optional',
  'editor',
  'nullable',
  'separators',
  'choices',
  'format',
  'min',
  'max',
  'caption',
  'description',
  'initialValue',
  'viewer',
  'units',
  'type',
  'optionsType',
  'step',
  'meta.url',
  'metaUrl',
];

const nonquotedValues = ['true', 'false'];

function buildStringOfOptions(input: any) {
  const optionsInString: string[] = [];
  const opt = input.options ?? {};
  let defaultValue = '';
  if (opt['initialValue'] && /[A-Za-z]/.test(opt['initialValue']) && !opt['initialValue'].includes('\'') && 
    !opt['initialValue'].includes('"') &&
    !nonquotedValues.includes(opt['initialValue'])!)
    opt['initialValue'] = `'${opt['initialValue']}'`;

  if (opt['initialValue']) 
    defaultValue = `= ${opt['initialValue']}`; 

  for (const [key, value] of Object.entries(opt)) {
    if (key === 'initialValue')
      continue;
    let val = value;
    let option = key;
    option = decoratorOptionToAnnotation.get(option) ?? option;

    if (Array.isArray(value))
      val = JSON.stringify(value);
    optionsInString.push(`${option}: ${val}`);
  }
  const optString = optionsInString.length> 0 ? `{ ${optionsInString.join('; ')} }`: '';
  return defaultValue? `${defaultValue} ${optString}` : `${optString}`;
}


interface ReservedDecorator{ 
  [decorator: string]: { metadata: FuncMetadata, genFunc: Function } 
}

export const reservedDecorators : ReservedDecorator = {
  viewer: {
    metadata: {
      role: FUNC_TYPES.VIEWER,
      inputs: [],
      outputs: [{name: 'result', type: 'viewer'}],
    },
    genFunc: generateClassFunc,
  },
  filter: {
    metadata: {
      role: FUNC_TYPES.FILTER,
      inputs: [],
      outputs: [{name: 'result', type: 'filter'}],
    },
    genFunc: generateClassFunc,
  },
  cellRenderer: {
    metadata: {
      role: FUNC_TYPES.CELL_RENDERER,
      inputs: [],
      outputs: [{name: 'renderer', type: 'grid_cell_renderer'}],
    },
    genFunc: generateClassFunc,
  },
  appTreeBrowser: {
    metadata: {
      tags: [],
      role: FUNC_TYPES.APP_TREE_BROWSER,
      inputs: [{type: 'dynamic', name: 'treeNode'}, {type: 'view', name: 'browseView'}],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  fileExporter: {
    metadata: {
      role: FUNC_TYPES.FILE_EXPORTER,
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  fileHandler: {
    metadata: {
      role: FUNC_TYPES.FILE_HANDLER,
      inputs: [{name: 'content', type: 'string'}],
      outputs: [{name: 'tables', type: 'list'}],
    },
    genFunc: generateFunc,
  },
  fileViewer: {
    metadata: {
      role: FUNC_TYPES.FILE_VIEWER,
      inputs: [{name: 'f', type: 'file'}],
      outputs: [{name: 'v', type: 'view'}],
    },
    genFunc: generateFunc,
  },
  settingsEditor: {
    metadata: {
      role: FUNC_TYPES.SETTINGS_EDITOR,
      inputs: [],
      outputs: [{name: 'result', type: 'widget'}],
    },
    genFunc: generateFunc,
  },
  func: {
    metadata: {
      tags: [],
      inputs: [],
      outputs: [{name: 'result', type: 'dynamic'}],
    },
    genFunc: generateFunc,
  },
  app: {
    metadata: {
      role: FUNC_TYPES.APP,
      inputs: [],
      outputs: [{name: 'result', type: 'view'}],
    },
    genFunc: generateFunc,
  },
  autostart: {
    metadata: {
      role: FUNC_TYPES.AUTOSTART,
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  init: {
    metadata: {
      role: FUNC_TYPES.INIT,
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  editor: {
    metadata: {
      role: FUNC_TYPES.EDITOR,
      inputs: [{name: 'call', type: 'funccall'}],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  panel: {
    metadata: {
      role: FUNC_TYPES.PANEL,
      inputs: [],
      outputs: [{name: 'result', type: 'widget'}],
    },
    genFunc: generateFunc,
  },
  folderViewer: {
    metadata: {
      role: FUNC_TYPES.FOLDER_VIEWER,
      inputs: [{name: 'folder', type: 'file'}, {name: 'files', type: 'list<file>'}],
      outputs: [{name: 'result', type: 'widget'}],
    },
    genFunc: generateFunc,
  },
  semTypeDetector: {
    metadata: {
      role: FUNC_TYPES.SEM_TYPE_DETECTOR,
      inputs: [{name: 'col', type: 'column'}],
      outputs: [{name: 'result', type: 'string'}],
    },
    genFunc: generateFunc,
  },
  dashboard: {
    metadata: {
      role: FUNC_TYPES.DASHBOARD,
      inputs: [],
      outputs: [{name: 'result', type: 'widget'}],
    },
    genFunc: generateFunc,
  },
  functionAnalysis: {
    metadata: {
      role: FUNC_TYPES.FUNCTION_ANALYSIS,
      inputs: [],
      outputs: [{name: 'result', type: 'view'}],
    },
    genFunc: generateFunc,
  },
  converter: {
    metadata: {
      role: FUNC_TYPES.CONVERTER,
      inputs: [{name: 'value', type: 'dynamic'}],
      outputs: [{name: 'result', type: 'dynamic'}],
    },
    genFunc: generateFunc,
  },
  demo: {
    metadata: {
      tags: [],
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  treeBrowser: {
    metadata: {
      tags: [],
      inputs: [{name: 'treeNode', type: 'dynamic'}, {name: 'browseView', type: 'view'}],
      outputs: [],
    },
    genFunc: generateFunc,
  },
  model: {
    metadata: {
      role: FUNC_TYPES.MODEL,
      inputs: [],
      outputs: [],
    },
    genFunc: generateFunc,
  },

};

/** Generates a DG function that instantiates a class. */
export function generateClassFunc(annotation: string, className: string, sep: string = '\n'): string {
  return annotation + `export function _${className}() {${sep}  return new ${className}();${sep}}${sep.repeat(2)}`;
}

export const primitives = new Set([
  'string',
  'string[]',
  'number',
  'number[]',
  'boolean',
  'boolean[]',
  'any',
  'Uint8Array',
  'Uint8Array[]',
  'void',
]);

/** Generates a DG function. */
export function generateFunc(
  annotation: string,
  funcName: string,
  sep: string = '\n',
  className: string = '',
  inputs: FuncParam[] = [],
  output: string,
  isAsync: boolean = false): string {
  // FuncCall can have an optional followed by mandatory args for ui reasons, typescript cannot
  let firstTsValidOptionalIdx = inputs.length-1;
  for (; firstTsValidOptionalIdx >= 0; firstTsValidOptionalIdx--) {
    const e = inputs[firstTsValidOptionalIdx];
    if (!e.optional)
      break;
  }
  const funcSigNature = (inputs.map((e, idx) => {
    const namePart = `${e.name}${(e.optional && idx >= firstTsValidOptionalIdx) ? '?': ''}`;
    // eslint-disable-next-line max-len
    let typePart = `${primitives.has(e.type ?? '') && !typesToAny.includes(e.type ?? '') ? e.type : (typesToAnnotation[e.type?.replace('[]', '') ?? ''] && !typesToAny.includes(e.type ?? '') ? e.type : 'any')}`;
    if (e.optional && idx < firstTsValidOptionalIdx)
      typePart += ' | undefined';
    return `${namePart}: ${typePart}`;
  })).join(', ');
  const funcArguments = (inputs.map((e) => e.name)).join(', ');

  const returnType = output ? ( primitives.has(output) ?
    (!isAsync? `: ${output} `: `: Promise<${output}> `) : (!isAsync? `: any `: `: Promise<any> `)) : '';
  // eslint-disable-next-line max-len
  return sep + annotation + `export ${isAsync ? 'async ' : ''}function ${funcName}(${funcSigNature}) ${returnType}{${sep}  ${output !== 'void'? 'return ': ''}${isAsync? 'await ': ''}${className.length > 0 ? `${className}.` : ''}${funcName}(${funcArguments});${sep}}${sep}`;
}

export function generateImport(className: string, path: string, sep: string = '\n'): string {
  return `import {${className}} from '${path}';${sep}`;
}

export function generateExport(className: string, sep: string = '\n'): string {
  return `export {${className}};${sep}`;
}
