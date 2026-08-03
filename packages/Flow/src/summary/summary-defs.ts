/** Curated, hand-written plain-language summaries for the functions a scientist
 *  actually chains, plus Flow's built-in node types — the data behind the
 *  heuristic node/flow summaries (U12 "the flow documents itself").
 *
 *  Coverage was chosen empirically from the live catalog (800 registered funcs:
 *  287 core, 122 Chem, 48 Bio, …) — see the retired summary-diag suite. Anything
 *  not listed here falls back to a humanized name in `summary-generator.ts`, so
 *  this map only needs the high-traffic functions to lift the average.
 *
 *  Keys are the **bare** DG function name, lower-cased (the generator strips any
 *  `Package:` prefix and lower-cases before lookup). Built-ins are keyed by their
 *  registered node type name. */

import {FlowNode} from '../rete/scheme';

export type SummaryTemplate = (node: FlowNode) => string;

// ---------- small formatters (shared with templates) ----------

export const str = (v: unknown): string => (v == null ? '' : String(v)).trim();
export const truncate = (s: string, n = 44): string => (s.length > n ? s.slice(0, n - 1) + '…' : s);
export const basename = (p: unknown): string => {
  const s = str(p);
  return s.split(/[\\/]/).filter(Boolean).pop() || s;
};
/** A node's literal input value (set in the panel; unset → undefined). */
export const iv = (node: FlowNode, key: string): unknown =>
  (node.inputValues as Record<string, unknown> | undefined)?.[key];
export const ivs = (node: FlowNode, key: string): string => str(iv(node, key));
/** A node property (paramName, outputType, …). */
export const prop = (node: FlowNode, key: string): string =>
  str((node.properties as Record<string, unknown> | undefined)?.[key]);

/** "MW, logP, HBA" from the boolean flags the user enabled on a properties node. */
function enabledFlags(node: FlowNode, flags: Array<[string, string]>): string {
  const on = flags.filter(([k]) => iv(node, k) === true || iv(node, k) === 'true').map(([, label]) => label);
  return on.join(', ');
}

// ---------- curated function summaries (bare name, lower-cased) ----------

export const CURATED_FUNC_SUMMARIES: Record<string, SummaryTemplate> = {
  // --- data sources (core) ---
  'openfile': (n) => {
    const f = basename(iv(n, 'fullPath'));
    return f ? `Loads file ${f}` : 'Loads a file';
  },
  'opentable': (n) => {
    const id = ivs(n, 'id');
    return id ? `Opens table ${id}` : 'Opens a table';
  },
  'dbquery': () => 'Runs a database query',
  'appendtables': () => 'Appends tables into one',
  'randomdata': (n) => `Generates random data${ivs(n, 'distribution') ? ` (${ivs(n, 'distribution')})` : ''}`,

  // --- combine ---
  'jointables': (n) => `Joins two tables${ivs(n, 'joinType') ? ` (${ivs(n, 'joinType')})` : ''}`,
  'linktables': () => 'Links two tables by key',
  'comparetables': () => 'Compares two tables',

  // --- transform / column ops (core) ---
  'addnewcolumn': (n) => {
    const name = ivs(n, 'name');
    const expr = ivs(n, 'expression');
    if (name && expr) return `Adds column “${name}” = ${expr}`;
    if (name) return `Adds column “${name}”`;
    return 'Adds a calculated column';
  },
  'addnewcolumnlist': () => 'Adds several calculated columns',
  'aggregate': () => 'Aggregates / groups the table',
  'unpivot': () => 'Unpivots columns into rows',
  'filterrows': () => 'Filters rows by a condition',
  'extractrows': () => 'Extracts rows matching a condition',
  'deleterows': () => 'Deletes rows matching a condition',
  'deletecolumns': () => 'Deletes columns',
  'keepcolumns': () => 'Keeps only selected columns',
  'extractcolumns': () => 'Extracts selected columns',
  'renamecolumn': (n) => {
    const a = ivs(n, 'oldColName'); const b = ivs(n, 'newColName');
    return a && b ? `Renames ${a} → ${b}` : 'Renames a column';
  },
  'changecolumnstype': (n) => `Changes column type${ivs(n, 'newType') ? ` → ${ivs(n, 'newType')}` : ''}`,
  'clonetable': () => 'Duplicates the table',
  'normalize': (n) => `Normalizes ${ivs(n, 'column') || 'a column'}${ivs(n, 'method') ? ` (${ivs(n, 'method')})` : ''}`,
  'missingvaluesimputation': () => 'Imputes missing values',
  'pca': (n) => `PCA${ivs(n, 'components') ? ` → ${ivs(n, 'components')} components` : ''}`,
  'cluster': (n) => `K-means clustering${ivs(n, 'clusters') ? ` (${ivs(n, 'clusters')})` : ''}`,
  'splitcolumn': () => 'Splits a column',

  // --- variables / plumbing (core; common in imported flows) ---
  'setvar': (n) => {
    const v = ivs(n, 'variableName') || str(n.label).replace(/^set:\s*/i, '');
    return v ? `Stores result as “${v}”` : 'Stores a variable';
  },
  'getvar': (n) => {
    const v = ivs(n, 'variableName') || str(n.label).replace(/^get:\s*/i, '');
    return v ? `Reads variable “${v}”` : 'Reads a variable';
  },

  // --- compute (core) ---
  'tostring': () => 'Converts a value to text',
  'concat': () => 'Concatenates values',

  // --- visualize (core) ---
  'createviewer': (n) => `Adds a ${ivs(n, 'viewerType') || ''} viewer`.replace('  ', ' '),
  'histogram': () => 'Adds a histogram',
  'scatterplot': () => 'Adds a scatter plot',
  'barchart': () => 'Adds a bar chart',
  'linechart': () => 'Adds a line chart',
  'piechart': () => 'Adds a pie chart',
  'boxplot': () => 'Adds a box plot',
  'treemap': () => 'Adds a tree map',
  'trellis': () => 'Adds a trellis plot',

  // --- Chem ---
  'addchempropertiescolumns': (n) => {
    const sel = enabledFlags(n, [
      ['MW', 'MW'], ['HBA', 'HBA'], ['HBD', 'HBD'], ['logP', 'logP'], ['logS', 'logS'],
      ['PSA', 'PSA'], ['rotatableBonds', 'rot. bonds'], ['stereoCenters', 'stereocenters'],
      ['moleculeCharge', 'charge'],
    ]);
    return sel ? `Computes ${sel}` : 'Computes chemical properties';
  },
  'calculatelogp': () => 'Calculates logP',
  'calculatelogd': () => 'Calculates logD',
  'calculatelogs': () => 'Calculates logS (solubility)',
  'calculatepka': () => 'Calculates pKa',
  'calculatepi': () => 'Calculates pI',
  'curate': () => 'Standardizes (curates) structures',
  'elementalanalysis': () => 'Elemental analysis of molecules',
  'murckoscaffolds': () => 'Extracts Murcko scaffolds',
  'generatescaffoldtree': () => 'Builds a scaffold tree',
  'filterbycatalogs': (n) => `Flags undesirable molecules${ivs(n, 'catalog') ? ` (${ivs(n, 'catalog')})` : ''}`,
  'findsimilar': () => 'Finds molecules similar to a query',
  'findmcs': () => 'Finds the maximum common substructure',
  'findrgroups': () => 'R-group decomposition',
  'findrgroupswithcore': () => 'R-group decomposition (given a core)',
  'chemspacetopmenu': () => 'Maps chemical space in 2D',
  'activitycliffs': () => 'Detects activity cliffs',
  'butinamoleculesclustering': () => 'Butina clustering of molecules',
  'bitbirchclusteringtopmenu': () => 'BitBIRCH clustering of molecules',
  'clustermcstopmenu': () => 'Common substructure per cluster',
  'generateconformers': () => 'Generates 3D conformers',
  'convertmolnotation': () => 'Converts molecule notation',
  'mutate': () => 'Generates analog structures',
  'deprotect': () => 'Removes protecting groups',

  // --- Bio ---
  'alignsequences': () => 'Multiple sequence alignment',
  'sequencespacetopmenu': () => 'Maps sequence space in 2D',
  'compositionanalysis': () => 'Sequence composition (WebLogo)',
  'getregion': (n) => {
    const a = ivs(n, 'start'); const b = ivs(n, 'end');
    return a || b ? `Extracts sequence region ${a || '?'}–${b || '?'}` : 'Extracts a sequence region';
  },
  'toatomiclevel': () => 'Converts sequences to atomic level',
  'sequenceidentityscoring': () => 'Scores sequence identity',
  'sequencesimilarityscoring': () => 'Scores sequence similarity',
  'seqidentity': () => 'Sequence identity to a reference',
  'similaritysearchtopmenu': () => 'Finds similar sequences',
  'diversitysearchtopmenu': () => 'Finds the most diverse sequences',
  'moleculestohelmtopmenu': () => 'Converts molecules to HELM',
  'splittomonomerstopmenu': () => 'Splits sequences into monomers',
  'importfasta': () => 'Opens a FASTA file',
};

// ---------- built-in node summaries (keyed by registered type name) ----------

export const BUILTIN_SUMMARIES: Record<string, SummaryTemplate> = {
  'Utilities/Select Column': (n) => `Picks column ${prop(n, 'columnName') || ivs(n, 'columnName') || ''}`.trim(),
  'Utilities/Select Columns': () => 'Picks several columns',
  'Utilities/Select Table': (n) => `Selects table ${prop(n, 'tableName') || str(n.label).replace(/^table:\s*/i, '') || ''}`.trim(),
  'Utilities/Add Table View': () => 'Opens the table in a grid',
  'Utilities/Log': () => 'Logs a message',
  'Utilities/Info': () => 'Shows an info message',
  'Utilities/Warning': () => 'Shows a warning',
  'Utilities/ToString': () => 'Converts a value to text',
  'Utilities/FromJSON': () => 'Parses JSON',
  'Utilities/ToJSON': () => 'Serializes to JSON',
  'Debug/Breakpoint': () => 'Pauses here when debugging',
};
