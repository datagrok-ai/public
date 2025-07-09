/* eslint-disable max-len */
/* eslint-disable no-tabs */
export const MIXFILE_VERSION = 1.00; // version number to use for newly created instances

export type MixfileMetadatum = string | (string | number)[];

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
	contents?:MixfileComponent[];
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
