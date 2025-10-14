import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {renderMolecule} from '../rendering/render-molecule';
/* eslint-disable max-len */
/* eslint-disable no-tabs */
export const MIXFILE_VERSION = 1.00; // version number to use for newly created instances

export type MixfileMetadatum = string | (string | number)[];

export const STRUCTURE_FIELDS = ['formula', 'molfile', 'inchi', 'inchiKey', 'smiles'];

export interface Mixfile extends MixfileComponent
{
	mixfileVersion:number;
}

export interface MixfileComponent
{
	name?:string;
	description?:string;
	synonyms?:string[];

	// molecular specification: none of them are mandatory, but if more than one is specified, they must refer to the same species,
	// to the extent that the format allows; InChI strings are expected to be standard, while SMILES are not
	formula?:string;
	molfile?:string;
	inchi?:string;
	inchiKey?:string;
	smiles?:string;

	// if the concentration is known, then these fields should be filled out as appropriate; if the concentration is a ratio,
	// it is relative to all of the components within the same branch
	quantity?:number | number[]; // a concentration numeric which is associated with the units below (two numbers in case of a range)
	error?:number; // optional standard error (applies to quantity when it's a scalar)
	ratio?:number[]; // a ratio, specified as [numerator, denominator]
	units?:string; // units for quantity (e.g. %, mol/L, g, etc.)
	relation?:string; // optional modifier when applied to quantity (e.g. >, <, ~)

	// identifiers that map the substance to external databases (e.g. PubChem, ChemSpider, CAS, etc.); identifiers are ID numbers, and
	// the meaning is implied by the context; links should be resolvable URLs, which are an alternative way of locating external
	// resources; in some cases where there are multiple identifiers, the value should be specified as an array
	identifiers?:Record<string, string | string[]>;
	links?:Record<string, string | string[]>;

	// metadata starts with an IRI that defines the core concept, and then may be followed by scalar data and/or other IRIs; this allows
	// discrete facts to be asserted about the component, as well as numeric values such as physical properties and their units
	metadata?:MixfileMetadatum[];

	// subcomponents: if this is a discrete molecular entity, then there will be none; usually there are either 0 or 2-or-more; in cases
	// where there are any subcomponents, any of the properties above apply to all of these subcomponents collectively
	contents?: MixfileComponent[];
}

// useful for cleaning up external JSON content
export const MIXFILE_COMPONENT_FIELDS =
[
  'name', 'description', 'synonyms', 'formula', 'molfile', 'inchi', 'inchiKey', 'smiles',
  'ratio', 'quantity', 'units', 'relation', 'identifiers', 'links', 'contents',
];
export const MIXFILE_ROOT_FIELDS =
[
  'mixfileVersion', ...MIXFILE_COMPONENT_FIELDS,
];

export interface MixtureWidgetObj {
  structures: string[];
  names: string[];
  quantities: string[];
  units: string[];
  ratios: string[];
  relations: string[];
  errors: (number | undefined)[];
  levels: number[];
  mixturesIds: string[];
}


export function transformMixfile(mixfile: Mixfile): MixtureWidgetObj {
  const structures: string[] = [];
  const names: string[] = [];
  const quantities: string[] = [];
  const units: string[] = [];
  const ratios: string[] = [];
  const relations: string[] = [];
  const errors: (number | undefined)[] = [];
  const levels: number[] = [];
  const mixturesIds: string[] = [];

  function traverse(component: MixfileComponent, depth: number, parentMixtureName: string = '') {
    const isPureContainer = !component.name && !component.molfile && !component.smiles &&
      !component.quantity && component.contents;

    if (!isPureContainer) {
      structures.push(component.molfile || component.smiles || '');
      names.push(component.name || '');
      quantities.push(
        component.quantity ?
          Array.isArray(component.quantity) ?
            `${component.quantity[0]}${component.quantity[1] ? ` - ${component.quantity[1]}` : ''}` :
            component.quantity.toString() : '',
      );
      units.push(component.units || '');
      ratios.push(
        component.ratio && component.ratio.length >= 2 ?
          `${component.ratio[0]}/${component.ratio[1]}` : '',
      );
      relations.push(component.relation || '');
      errors.push(component.error);
      levels.push(depth);
      mixturesIds.push(parentMixtureName);
    }

    // If this node is a mixture (has contents and a name), pass its name as parent for its children
    const nextParentMixtureName = (component.contents && component.name) ? component.name : '';
    if (component.contents)
      component.contents.forEach((child) => traverse(child, depth + 1, nextParentMixtureName));
  }

  if (Array.isArray(mixfile.contents)) {
    const rootName = mixfile.name || '';
    mixfile.contents.forEach((child) => traverse(child, 1, rootName));
  }

  return {structures, names, quantities, units, ratios, relations, errors, levels, mixturesIds};
}

export async function createMixtureWidget(mixture: string): Promise<DG.Widget> {
  const mixtureObj = JSON.parse(mixture) as Mixfile;

  const mw = transformMixfile(mixtureObj);
  let df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('structure', mw.structures),
    DG.Column.fromStrings('name', mw.names),
    DG.Column.fromStrings('relation', mw.relations),
    DG.Column.fromStrings('quantity', mw.quantities),
    DG.Column.fromStrings('units', mw.units),
    DG.Column.fromStrings('ratio', mw.ratios),
    DG.Column.fromList(DG.TYPE.FLOAT, 'SE', mw.errors),
    DG.Column.fromList(DG.TYPE.INT, 'level', mw.levels),
    DG.Column.fromStrings('parent mixture name', mw.mixturesIds),
  ]);
  if (!df)
    df = DG.DataFrame.create();
  else
    await grok.data.detectSemanticTypes(df);
  const grid = df.plot.grid();
  grid.root.style.height = '300px';

  const addToWorkspaceButton = ui.icons.add(async () => {
    grok.shell.addTableView(df);
  }, 'Add to workspace');
  addToWorkspaceButton.classList.add('chem-mixture-widget-add-table-to-workspace');

  return new DG.Widget(ui.divV([
    addToWorkspaceButton,
    grid.root,
  ]));
}


export function createComponentPane(component: MixfileComponent): HTMLElement {
  const fieldsForTableFromMap: {[key: string]: any} = {};
  const keys = Object.keys(component);
  const accordions = ui.divV([]);
  const structure = component.molfile ?? component.smiles ?? '';
  for (const key of keys) {
    //exclude structure fields
    if (STRUCTURE_FIELDS.includes(key))
      continue;
    //handling metadata separately, since it is an array of either scalar or array types
    if (key === 'metadata') {
      let str = '';
      for (const val of (component as any)[key])
        str += Array.isArray(val) ? val.join(',') : val;

      fieldsForTableFromMap[key] = str;
    } else if (Array.isArray((component as any)[key]) && (component as any)[key].length && //for array of scalar types () - join with comma
        (typeof (component as any)[key][0] === 'string' || typeof (component as any)[key][0] === 'number'))
      fieldsForTableFromMap[key] = (component as any)[key].join(',');
    //dictionaries (identifiers, links)
    else if (!Array.isArray((component as any)[key]) && typeof (component as any)[key] === 'object' && (component as any)[key] != null) {
      const innerDictAcc = ui.accordion(key);
      const obj = (component as any)[key];
      innerDictAcc.addPane(key, () => {
        const fields: {[key: string]: any} = {};
        Object.keys(obj).forEach((k: string) => {
          if (Array.isArray(obj[k]))
            fields[k] = obj[k].join(',');
          else
            fields[k] = obj[k];
        });
        return ui.tableFromMap(fields);
      });
      accordions.append(innerDictAcc.root);
    } else if (key === 'contents') {
      const contentsAcc = ui.accordion(key);
      contentsAcc.addPane(key, () => {
        const innerContentsAcc = ui.accordion();
        for (let i = 0; i < (component as any)[key].length; i++)
          innerContentsAcc.addPane((component as any)[key][i].name ?? `component ${i + 1}`, () => createComponentPane((component as any)[key][i]));

        return innerContentsAcc.root;
      });
      accordions.append(contentsAcc.root);
    } else //scalar types
      fieldsForTableFromMap[key] = (component as any)[key];
  }
  //handling structure fileds separately (show just one structure field in case molfile and smiles are missing)
  addStructureFields(fieldsForTableFromMap, component);
  const resDiv = ui.divV([]);
  if (structure)
    resDiv.append(renderMolecule(structure, {renderer: 'RDKit'}));
  resDiv.append(ui.tableFromMap(fieldsForTableFromMap));
  resDiv.append(accordions);
  return resDiv;
}

export function addStructureFields(dict: {[key: string]: any}, comp: MixfileComponent) {
  //use only one of the structure fields in case molfile and smiles are missing
  if (!comp.molfile && !comp.smiles) {
    if (comp.inchi)
      dict['inchi'] = comp.inchi;
    else if (comp.inchiKey)
      dict['inchiKey'] = comp.inchiKey;
    else if (comp.formula)
      dict['formula'] = comp.formula;
  }
}
