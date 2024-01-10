/* eslint-disable no-unused-vars */
/* eslint-disable valid-jsdoc */
import {FuncMetadata} from './interfaces';

export const headerParams = ['name', 'description', 'tags', 'inputs', 'outputs'];

export enum pseudoParams {
  EXTENSION = 'extension',
  EXTENSIONS = 'extensions',
  INPUT_TYPE = 'inputType',
}

export enum FUNC_TYPES {
  CELL_RENDERER = 'cellRenderer',
  FILE_EXPORTER = 'fileExporter',
  FILE_IMPORTER = 'file-handler',
  FILE_VIEWER = 'fileViewer',
  SETTINGS_EDITOR = 'packageSettingsEditor',
  VIEWER = 'viewer',
  FILTER = 'filter',
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
      data.tags.concat(data[pseudoParams.EXTENSIONS].map((ext: string) => 'fileViewer-' + ext)).join() :
      data.tags.join()}${sep}`;
  }
  if (data.inputs) {
    for (const input of data.inputs) {
      s += comment + 'input: ' + (isFileImporter && data[pseudoParams.INPUT_TYPE] ?
        data[pseudoParams.INPUT_TYPE] : input.type) + (input.name ? ` ${input.name}` : '') + sep;
    }
  }
  if (data.outputs) {
    for (const output of data.outputs) 
      s += comment + 'output: ' + output.type + (output.name ? ` ${output.name}` : '') + sep;
    
  }
  for (const parameter in data) {
    if (parameter === pseudoParams.EXTENSION || parameter === pseudoParams.INPUT_TYPE)
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

export const reservedDecorators: {[decorator: string]: {metadata: FuncMetadata, genFunc: Function}} = {
  viewer: {
    metadata: {
      tags: [FUNC_TYPES.VIEWER],
      inputs: [],
      outputs: [{name: 'result', type: 'viewer'}],
    },
    genFunc: generateClassFunc,
  },
  filter: {
    metadata: {
      tags: [FUNC_TYPES.FILTER],
      inputs: [],
      outputs: [{name: 'result', type: 'filter'}],
    },
    genFunc: generateClassFunc,
  },
  cellRenderer: {
    metadata: {
      tags: [FUNC_TYPES.CELL_RENDERER],
      inputs: [],
      outputs: [{name: 'renderer', type: 'grid_cell_renderer'}],
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
      inputs: [{name: 'content', type: 'string'}],
      outputs: [{name: 'tables', type: 'list'}],
    },
    genFunc: generateFunc,
  },
  fileViewer: {
    metadata: {
      tags: [FUNC_TYPES.FILE_VIEWER],
      inputs: [{name: 'f', type: 'file'}],
      outputs: [{name: 'v', type: 'view'}],
    },
    genFunc: generateFunc,
  },
  settingsEditor: {
    metadata: {
      tags: [FUNC_TYPES.SETTINGS_EDITOR],
      inputs: [],
      outputs: [{name: 'result', type: 'widget'}],
    },
    genFunc: generateFunc,
  },
};

/** Generates a DG function that instantiates a class. */
export function generateClassFunc(annotation: string, className: string, sep: string = '\n'): string {
  return annotation + `export function _${className}() {${sep}  return new ${className}();${sep}}${sep.repeat(2)}`;
}

/** Generates a DG function. */
export function generateFunc(annotation: string, funcName: string, sep: string = '\n'): string {
  return annotation + `export function _${funcName}() {${sep}  return ${funcName}();${sep}}${sep.repeat(2)}`;
}

export function generateImport(className: string, path: string, sep: string = '\n'): string {
  return `import {${className}} from '${path}';${sep}`;
}

export function generateExport(className: string, sep: string = '\n'): string {
  return `export {${className}};${sep}`;
}
