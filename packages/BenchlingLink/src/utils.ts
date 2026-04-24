import * as DG from 'datagrok-api/dg';

// Translate TS-friendly underscore suffixes (_anyOf, _anyOf_caseSensitive, _gt/_gte/_lt/_lte,
// _mol, _smiles) into Benchling's dotted query-param names (names.anyOf, emptyPositions.gt, …).
// Only the documented suffixes are mapped; ordinary keys pass through unchanged.
const SUFFIX_TRANSLATIONS: Array<[RegExp, string]> = [
  [/_anyOf_caseSensitive$/, '.anyOf.caseSensitive'],
  [/_anyOf$/, '.anyOf'],
  [/_gte$/, '.gte'],
  [/_gt$/, '.gt'],
  [/_lte$/, '.lte'],
  [/_lt$/, '.lt'],
  [/_mol$/, '.mol'],
  [/_smiles$/, '.smiles'],
];

function translateParamName(key: string): string {
  for (const [re, repl] of SUFFIX_TRANSLATIONS) {
    if (re.test(key))
      return key.replace(re, repl);
  }
  return key;
}

export function buildQueryString(params: Record<string, any>): string {
  const esc = encodeURIComponent;
  return Object.entries(params)
    .filter(([_, v]) => {
      if (v === undefined || v === null || v === '') return false;
      if (Array.isArray(v) && v.length === 0) return false;
      return true;
    })
    .map(([k, v]) => `${esc(translateParamName(k))}=${esc(Array.isArray(v) ? v.join(',') : v)}`)
    .join('&');
}

export function randomDnaSequence(length: number): string {
  const dna = 'ATGC';
  let seq = '';
  for (let i = 0; i < length; i++)
    seq += dna[Math.floor(Math.random() * dna.length)];
  return seq;
}

function inferColumnType(value: any): DG.COLUMN_TYPE {
  if (typeof value === 'boolean') return DG.COLUMN_TYPE.BOOL;
  if (typeof value === 'bigint') return DG.COLUMN_TYPE.BIG_INT;
  if (typeof value === 'number')
    return Number.isInteger(value) ? DG.COLUMN_TYPE.INT : DG.COLUMN_TYPE.FLOAT;
  return DG.COLUMN_TYPE.STRING;
}

export function dataFrameFromObjects(objects: Record<string, any>[]): DG.DataFrame {
  if (!objects || objects.length === 0)
    return DG.DataFrame.create();

  // Helper function to process array of objects - extract id or identifier filed and concatenate
  const processObjectArray = (arr: any[]): string => {
    if (!Array.isArray(arr)) return '';
    
    return arr.map(obj => {
      if (typeof obj !== 'object' || obj === null) return String(obj);
      
      // Try to find identifier field
      let identifier = obj.id;
      if (identifier === undefined) {
        const identifierKey = Object.keys(obj).find(key => 
          key.toLowerCase().includes('identifier')
        );
        identifier = identifierKey ? obj[identifierKey] : undefined;
      }
      
      return identifier !== undefined ? String(identifier) : '';
    }).filter(Boolean).join(',');
  };

  // Helper function to flatten objects recursively
  const flattenObject = (obj: Record<string, any>, prefix = ''): Record<string, any> => {
    const result: Record<string, any> = {};
    
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        const value = obj[key];
        const newKey = prefix ? `${prefix}.${key}` : key;
        
        if (Array.isArray(value) && value.some((v) => v !== null && typeof v === 'object')) {
          // Handle array of objects specially
          result[newKey] = processObjectArray(value);
        }
        else if (typeof value === 'object' && value !== null && !Array.isArray(value)) {
          // Recursively flatten nested objects
          Object.assign(result, flattenObject(value, newKey));
        } else {
          // Convert regular arrays to comma-separated strings
          result[newKey] = Array.isArray(value) ? value.join(',') : value;
        }
      }
    }
    
    return result;
  };

  // First flatten all objects in the array
  const flattenedObjects = objects.map(obj => flattenObject(obj));

  // Now get all unique keys from all flattened objects
  const allKeys = new Set<string>();
  for (const obj of flattenedObjects) {
    Object.keys(obj).forEach(key => allKeys.add(key));
  }

  // Create columns for each key
  const columns: DG.Column[] = [];
  for (const key of Array.from(allKeys)) {
    // Find first non-null, non-undefined value to infer type
    let firstValue = flattenedObjects.find(obj => obj[key] !== undefined && obj[key] !== null)?.[key];
    if (firstValue === undefined || firstValue === null) firstValue = '';

    const colType = inferColumnType(firstValue);
    columns.push(DG.Column.fromList(colType, key,
      flattenedObjects.map(obj => obj[key] !== undefined && obj[key] !== null ? obj[key] : null)));
  }

  return DG.DataFrame.fromColumns(columns);
}