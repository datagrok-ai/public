export const headerParams = ['name', 'description', 'tags', 'inputs', 'outputs'];

/** Generates an annotation header for a function based on provided metadata. */
export function getFuncAnnotation(data: FuncMetadata, comment: string = '//', sep: string = '\n'): string {
  let s = '';
  if (data.name)
    s += `${comment}name: ${data.name}${sep}`;
  if (data.description)
    s += `${comment}description: ${data.description}${sep}`;
  if (data.tags)
    s += `${comment}tags: ${data.tags.join()}${sep}`;
  if (data.inputs) {
    for (const input of data.inputs) {
      s += comment + 'input: ' + input.type + (input.name ? ` ${input.name}` : '') + sep;
    }
  }
  if (data.outputs) {
    for (const output of data.outputs) {
      s += comment + 'output: ' + output.type + (output.name ? ` ${output.name}` : '') + sep;
    }
  }
  for (const parameter in data) {
    if (!headerParams.includes(parameter)) {
      s += `${comment}meta.${parameter}: ${data[parameter]}${sep}`;
    }
  }
  return s;
}

export const reservedDecorators: {[decorator: string]: {metadata: FuncMetadata, genFunc: Function}} = {
  grokViewer: {
    metadata: {
      tags: ['viewer'],
      inputs: [],
      outputs: [{name: 'v', type: 'viewer'}],
    },
    genFunc: generateViewerFunc,
  },
  // TODO: add other decorators
};

/** Generates a DG function. */
export function generateViewerFunc(annotation: string, className: string, sep: string = '\n'): string {
  // TODO: add an import statement for the class
  return annotation + `export function _${className}() {${sep}  return new ${className}();${sep}}${sep.repeat(2)}`;
}

export function generateImport(className: string, path: string, sep: string = '\n'): string {
  return `import {${className}} from '${path}';${sep}`;
}

export interface Indexable {
  [key: string]: any,
}

export interface FuncParam {
  name?: string,
  type: string,
}

export interface FuncMetadata extends Indexable {
  name?: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
}
