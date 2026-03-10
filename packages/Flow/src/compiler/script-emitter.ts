import {LGraph} from 'litegraph.js';
import {CompiledStep, compileGraph} from './graph-compiler';
import {FuncFlowNode} from '../types/funcflow-node';

export interface ScriptSettings {
  name: string;
  description: string;
  tags: string[];
}

/** Generates a valid Datagrok script from compiled steps */
export function emitScript(graph: LGraph, settings: ScriptSettings): string {
  const steps = compileGraph(graph);
  const lines: string[] = [];

  // --- Header ---
  lines.push(`//name: ${settings.name}`);
  if (settings.description)
    lines.push(`//description: ${settings.description}`);
  lines.push('//language: javascript');
  if (settings.tags.length > 0)
    lines.push(`//tags: ${settings.tags.join(', ')}`);

  // --- Input declarations (table/dataframe inputs first) ---
  const inputSteps = steps.filter((s) => s.nodeType === 'input');
  inputSteps.sort((a, b) => {
    const aNode = graph.getNodeById(a.nodeId) as FuncFlowNode | undefined;
    const bNode = graph.getNodeById(b.nodeId) as FuncFlowNode | undefined;
    const aIsTable = aNode?.dgOutputType === 'dataframe';
    const bIsTable = bNode?.dgOutputType === 'dataframe';
    if (aIsTable && !bIsTable) return -1;
    if (!aIsTable && bIsTable) return 1;
    return 0;
  });
  for (const step of inputSteps) {
    const node = graph.getNodeById(step.nodeId) as FuncFlowNode | undefined;
    if (!node) continue;
    const inputLine = buildInputLine(step, node);
    if (inputLine) lines.push(inputLine);
  }

  // --- Output declarations ---
  const outputSteps = steps.filter((s) => s.nodeType === 'output');
  for (const step of outputSteps) {
    const node = graph.getNodeById(step.nodeId) as FuncFlowNode | undefined;
    if (!node) continue;
    const outputLine = buildOutputLine(step, node);
    if (outputLine) lines.push(outputLine);
  }

  lines.push('');

  // --- Body ---
  for (const step of steps) {
    if (step.nodeType === 'input') continue;
    if (step.nodeType === 'output') {
      const inputExpr = step.inputs.get('value') || 'undefined';
      if (inputExpr !== 'undefined')
        lines.push(`${step.variableName} = ${inputExpr};`);
      continue;
    }

    if (step.nodeType === 'utility') {
      const codeLine = emitUtilityStep(step);
      if (codeLine) lines.push(codeLine);
      continue;
    }

    lines.push(emitFuncStep(step));
  }

  return lines.join('\n');
}

function buildInputLine(step: CompiledStep, node: FuncFlowNode): string | null {
  const paramName = step.properties['paramName'];
  const outputType = node.dgOutputType || 'dynamic';
  const description = step.properties['description'] || '';

  let line = `//input: ${outputType} ${paramName}`;

  const qualifiers: string[] = [];

  const defaultVal = step.properties['defaultValue'];
  if (defaultVal !== undefined && defaultVal !== '' && defaultVal !== null)
    line = `//input: ${outputType} ${paramName} = ${formatHeaderDefault(defaultVal, outputType)}`;

  if (step.properties['typeFilter'])
    qualifiers.push(`type: ${step.properties['typeFilter']}`);
  if (step.properties['semTypeFilter'])
    qualifiers.push(`semType: ${step.properties['semTypeFilter']}`);
  if (step.properties['semType'])
    qualifiers.push(`semType: ${step.properties['semType']}`);
  if (step.properties['nullable'] === true)
    qualifiers.push('nullable: true');
  if (step.properties['caption'])
    qualifiers.push(`caption: ${step.properties['caption']}`);
  if (step.properties['choices']) {
    const items = step.properties['choices'].split(',').map((s: string) => s.trim()).filter(Boolean);
    qualifiers.push(`choices: [${items.map((s: string) => `"${s}"`).join(', ')}]`);
  }
  if (step.properties['min'] !== undefined && step.properties['min'] !== '')
    qualifiers.push(`min: ${step.properties['min']}`);
  if (step.properties['max'] !== undefined && step.properties['max'] !== '')
    qualifiers.push(`max: ${step.properties['max']}`);
  if (step.properties['showSlider'] === true)
    qualifiers.push('showSlider: true');

  if (qualifiers.length > 0)
    line += ` {${qualifiers.join('; ')}}`;

  if (description)
    line += ` [${description}]`;

  return line;
}

function buildOutputLine(step: CompiledStep, node: FuncFlowNode): string | null {
  const paramName = step.properties['paramName'];
  // For Value Output nodes, use the user-selected outputType property;
  // for Table Output, use the fixed dgOutputType ('dataframe')
  const outputType = step.properties['outputType'] || node.dgOutputType || 'dynamic';
  return `//output: ${outputType} ${paramName}`;
}

function emitFuncStep(step: CompiledStep): string {
  const params: string[] = [];
  for (const [name, expr] of step.inputs.entries())
    params.push(`${name}: ${expr}`);
  const paramsStr = params.length > 0 ? `{${params.join(', ')}}` : '{}';
  return `let ${step.variableName} = await grok.functions.call('${step.funcName}', ${paramsStr});`;
}

function emitUtilityStep(step: CompiledStep): string | null {
  switch (step.funcName) {
  case 'Select Column': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    const colName = step.properties['columnName'] || '';
    return `let ${step.variableName} = ${tableExpr}.col('${colName}');`;
  }
  case 'Select Columns': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    const colNames = (step.properties['columnNames'] || '').split(',').map((n: string) => n.trim()).filter(Boolean);
    const colExprs = colNames.map((n: string) => `${tableExpr}.col('${n}')`).join(', ');
    return `let ${step.variableName} = [${colExprs}];`;
  }
  case 'Add Table View': {
    const tableExpr = step.inputs.get('table') || 'undefined';
    return `let ${step.variableName} = grok.shell.addTableView(${tableExpr});`;
  }
  case 'Log': {
    const valueExpr = step.inputs.get('value') || 'undefined';
    const label = step.properties['label'] || '';
    return label ?
      `console.log('${label}:', ${valueExpr});` :
      `console.log(${valueExpr});`;
  }
  case 'Info': {
    const msgExpr = step.inputs.get('message') || '\'\'';
    return `grok.shell.info(${msgExpr});`;
  }
  case 'Warning': {
    const msgExpr = step.inputs.get('message') || '\'\'';
    return `grok.shell.warning(${msgExpr});`;
  }
  case 'ToString': {
    const valueExpr = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = (${valueExpr}).toString();`;
  }
  case 'String': {
    const val = step.properties['value'] ?? '';
    return `let ${step.variableName} = ${JSON.stringify(String(val))};`;
  }
  case 'Int': {
    const val = Math.round(Number(step.properties['value'] ?? 0));
    return `let ${step.variableName} = ${val};`;
  }
  case 'Double': {
    const val = Number(step.properties['value'] ?? 0);
    return `let ${step.variableName} = ${val};`;
  }
  case 'Boolean': {
    const val = step.properties['value'] ? 'true' : 'false';
    return `let ${step.variableName} = ${val};`;
  }
  case 'List': {
    const raw = String(step.properties['value'] ?? '');
    const items = raw.split(',').map((s: string) => s.trim()).filter(Boolean);
    const arrStr = items.map((s: string) => JSON.stringify(s)).join(', ');
    return `let ${step.variableName} = [${arrStr}];`;
  }
  case 'Equals (==)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) == (${b});`;
  }
  case 'Not Equals (!=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) != (${b});`;
  }
  case 'Greater Than (>)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) > (${b});`;
  }
  case 'Greater Or Equal (>=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) >= (${b});`;
  }
  case 'Less Than (<)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) < (${b});`;
  }
  case 'Less Or Equal (<=)': {
    const a = step.inputs.get('a') || 'undefined';
    const b = step.inputs.get('b') || 'undefined';
    return `let ${step.variableName} = (${a}) <= (${b});`;
  }
  case 'Contains': {
    const text = step.inputs.get('text') || '\'\'';
    const substring = step.inputs.get('substring') || '\'\'';
    return `let ${step.variableName} = (${text}).includes(${substring});`;
  }
  case 'Starts With': {
    const text = step.inputs.get('text') || '\'\'';
    const prefix = step.inputs.get('prefix') || '\'\'';
    return `let ${step.variableName} = (${text}).startsWith(${prefix});`;
  }
  case 'Ends With': {
    const text = step.inputs.get('text') || '\'\'';
    const suffix = step.inputs.get('suffix') || '\'\'';
    return `let ${step.variableName} = (${text}).endsWith(${suffix});`;
  }
  case 'Is Null': {
    const val = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = (${val}) == null;`;
  }
  case 'FromJSON': {
    const jsonExpr = step.inputs.get('json') || '\'{}\'';
    return `let ${step.variableName} = JSON.parse(${jsonExpr});`;
  }
  case 'ToJSON': {
    const valExpr = step.inputs.get('value') || 'undefined';
    return `let ${step.variableName} = JSON.stringify(${valExpr});`;
  }
  default:
    return `// Unknown utility: ${step.funcName}`;
  }
}

function formatHeaderDefault(value: any, type: string): string {
  switch (type) {
  case 'string':
    return `"${String(value)}"`;
  case 'bool':
    return value ? 'true' : 'false';
  case 'int':
  case 'double':
  case 'num':
    return String(value);
  default:
    return String(value);
  }
}
