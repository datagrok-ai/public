/** Built-in node catalog — one shared definition so the toolbox sections and the
 *  drag-out suggestion menu search the same texts. */

export interface BuiltinNodeDef {
  name: string;
  type: string;
  desc: string;
}

export interface BuiltinSection {
  title: string;
  nodes: BuiltinNodeDef[];
  tip: string;
}

export const BUILTIN_SECTIONS: BuiltinSection[] = [
  {
    title: 'Inputs',
    tip: 'Script input parameters (become //input: lines)',
    nodes: [
      {name: 'Property Input', type: 'Inputs/Property Input',
        desc: 'Mimics the function input it connects to — type, editor, and qualifiers come from the target parameter'},
      {name: 'Table Input', type: 'Inputs/Table Input', desc: 'Dataframe input parameter'},
      {name: 'Column Input', type: 'Inputs/Column Input', desc: 'Single column from a table'},
      {name: 'Column List Input', type: 'Inputs/Column List Input', desc: 'Multiple columns from a table'},
      {name: 'String Input', type: 'Inputs/String Input', desc: 'Text input parameter'},
      {name: 'Sketcher Input', type: 'Inputs/Sketcher Input', desc: 'A molecule, sketched — a string parameter tagged semType: Molecule'},
      {name: 'Helm Input', type: 'Inputs/Helm Input', desc: 'A macromolecule, edited in the HELM editor — a string parameter tagged semType: Macromolecule'},
      {name: 'Number Input', type: 'Inputs/Number Input', desc: 'Floating-point number input'},
      {name: 'Int Input', type: 'Inputs/Int Input', desc: 'Integer number input'},
      {name: 'Boolean Input', type: 'Inputs/Boolean Input', desc: 'True/false toggle input'},
      {name: 'DateTime Input', type: 'Inputs/DateTime Input', desc: 'Date and time input'},
      {name: 'File Input', type: 'Inputs/File Input', desc: 'File upload input'},
      {name: 'Map Input', type: 'Inputs/Map Input', desc: 'Key-value map input'},
      {name: 'Dynamic Input', type: 'Inputs/Dynamic Input', desc: 'Dynamically-typed input'},
      {name: 'String List Input', type: 'Inputs/String List Input', desc: 'List of strings input'},
      {name: 'Blob Input', type: 'Inputs/Blob Input', desc: 'Binary data input'},
    ],
  },
  {
    title: 'Outputs',
    tip: 'Script output parameters (become //output: lines)',
    nodes: [
      {name: 'Table Output', type: 'Outputs/Table Output', desc: 'Marks a dataframe as script output'},
      {name: 'Value Output', type: 'Outputs/Value Output', desc: 'Marks a value as script output (configurable type)'},
    ],
  },
  {
    title: 'Constants',
    tip: 'Constant literal values',
    nodes: [
      {name: 'String', type: 'Constants/String', desc: 'A constant text value'},
      {name: 'Int', type: 'Constants/Int', desc: 'A constant integer value'},
      {name: 'Double', type: 'Constants/Double', desc: 'A constant floating-point value'},
      {name: 'Boolean', type: 'Constants/Boolean', desc: 'A constant true/false value'},
      {name: 'List', type: 'Constants/List', desc: 'A constant list of comma-separated values'},
    ],
  },
  {
    title: 'Utilities',
    tip: 'Helper operations (logging, type conversion, etc.)',
    nodes: [
      {name: 'Select Column', type: 'Utilities/Select Column', desc: 'Gets a column from a table by name'},
      {name: 'Select Columns', type: 'Utilities/Select Columns', desc: 'Gets multiple columns by names (comma-separated)'},
      {name: 'Select Table', type: 'Utilities/Select Table', desc: 'Gets an open table by name via grok.shell.tableByName()'},
      {name: 'Add Table View', type: 'Utilities/Add Table View', desc: 'Opens a table in a new view via grok.shell.addTableView()'},
      {name: 'Log', type: 'Utilities/Log', desc: 'Logs a value to the browser console'},
      {name: 'Info', type: 'Utilities/Info', desc: 'Shows an info balloon via grok.shell.info()'},
      {name: 'Warning', type: 'Utilities/Warning', desc: 'Shows a warning balloon via grok.shell.warning()'},
      {name: 'ToString', type: 'Utilities/ToString', desc: 'Converts any value to a string'},
      {name: 'To Semantic Value', type: 'Utilities/To Semantic Value',
        desc: 'Wraps a value into a semantic value (DG.SemanticValue) of the given semantic type, e.g. Molecule'},
      {name: 'FromJSON', type: 'Utilities/FromJSON', desc: 'Parses a JSON string into an object'},
      {name: 'ToJSON', type: 'Utilities/ToJSON', desc: 'Serializes a value to a JSON string'},
    ],
  },
  {
    title: 'Debug',
    tip: 'Debugging and execution control nodes',
    nodes: [
      {name: 'Breakpoint', type: 'Debug/Breakpoint', desc: 'Pauses execution in debug mode until Continue is clicked'},
    ],
  },
];

let _descByType: Map<string, string> | null = null;

/** Description of a built-in node type, '' for unknown / DG-function types. */
export function builtinNodeDesc(typeName: string): string {
  if (!_descByType) {
    _descByType = new Map();
    for (const s of BUILTIN_SECTIONS)
      for (const n of s.nodes) _descByType.set(n.type, n.desc);
  }
  return _descByType.get(typeName) ?? '';
}
