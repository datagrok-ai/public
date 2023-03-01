/** Usage examples:
 * 
 * @grokViewer()
 * class TestViewer {
 *   constructor() {
 *     console.log('Viewer constructed');
 *   }
 * }
 * 
 * @grokViewer({
 *   name: 'Test Viewer',
 *   description: 'Creates a Test Viewer instance',
 *   icon: 'images/icon.png',
 *   toolbox: true,
 * })
 * class TestViewer {
 *   constructor() {
 *     console.log('Viewer constructed');
 *   }
 * }
 */
export function grokViewer(options?: {
  name?: string,
  description?: string,
  icon?: string,
  toolbox?: boolean,
}) {
  const funcMetadata = {
    tags: ['viewer'],
    inputs: [],
    outputs: [{name: '', type: 'viewer'}],
    ...options,
  };
  const annotation = getFuncAnnotation(funcMetadata);
  return function(constructor: Function) {
    // TODO: export the generated function to package.ts
    console.log(annotation +
      `export function _${constructor.name}() {\n  return new ${constructor.name}();\n}`);
  };
}

function getFuncAnnotation(data: FuncMetadata, comment: string = '//', sep: string = '\n'): string {
  let s = '';
  const parameters = ['name', 'description', 'tags', 'inputs', 'outputs'];
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
    if (!parameters.includes(parameter)) {
      s += `${comment}meta.${parameter}: ${data[parameter]}${sep}`;
    }
  }

  return s;
}

interface Indexable {
  [key: string]: any,
}

interface FuncMetadata extends Indexable {
  name?: string,
  description?: string,
  tags?: string[],
  inputs: FuncParam[],
  outputs: FuncParam[], 
}

interface FuncParam {
  name?: string,
  type: string,
}
