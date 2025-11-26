/* eslint-disable camelcase */
/* eslint-disable max-len */
import {Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {_package} from '../../package';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {HELM_RGROUP_FIELDS as RGP} from '@datagrok-libraries/bio/src/utils/const';
import * as grok from 'datagrok-api/grok';
import {getCorrectedMolBlock} from './monomer-manager/monomer-manager';

/**
 *  Exaple r groups
 * {
        "capGroupSMILES": "[*:1][H]",
        "alternateId": "R1-H",
        "capGroupName": "H",
        "label": "R1"
      },
      {
        "capGroupSMILES": "O[*:2]",
        "alternateId": "R2-OH",
        "capGroupName": "OH",
        "label": "R2"
      }

 */

const cx_smiles_regexp = /.*\|\$.*R.*\$\|/;
const rgroup_regexp = /\[R(\d+)\]/;
const rgroup_regexpg = /\[R(\d+)\]/g;
const ambig_regexp = /\[\*:(\d+)\]/g;

export type MonomerWithoutSymbol = Omit<Monomer, 'symbol'>;

export function getMonomerFromRSmiles(smiles: string, polymerType?: PolymerType): MonomerWithoutSymbol | null {
  const rgroupNumbers = Array.from(smiles.matchAll(rgroup_regexpg)).map((m) => m[1]);
  const res: MonomerWithoutSymbol = {
    name: 'Explicit SMILES Monomer',
    smiles: smiles,
    polymerType: polymerType ?? 'CHEM',
    molfile: '',
    rgroups: rgroupNumbers.map((numString) => ({
      [RGP.LABEL]: `R${numString}`,
      [RGP.CAP_GROUP_NAME]: `H`,
      [RGP.CAP_GROUP_SMILES]: `[*:${numString}][H]`,
      [RGP.ALTERNATE_ID]: `R${numString}-H`,
    })),
    author: 'Datagrok auto-generated',
    id: 0,
    createDate: null,
    monomerType: 'Backbone',
  };

  try {
    //try to generate corrected molfile and smiles
    let corSmiles = smiles;
    res.rgroups.forEach((rg) => {
      const labelNum = rg[RGP.LABEL].substring(1); // R1 -> 1
      corSmiles = corSmiles.replace(`[R${labelNum}]`, `[*:${labelNum}]`);
    });
    const molFile = getCorrectedMolBlock(grok.chem.convert(corSmiles, grok.chem.Notation.Smiles, grok.chem.Notation.MolBlock));
    res.molfile = molFile;
    res.smiles = corSmiles;
  } catch (e) {
    _package.logger.error(`getMonomerFromRSmiles: cannot convert SMILES to Molfile: ${smiles}\n${e}`);
  }

  return res;
}

/** Generate Monomer Object directly from inline smiles
 * Purely string based, no external calls
 *
 * Currently accepts (to be extended):
 *
 * cxsmiles written as *N[C@H](C(=O)*)Cc1ccc(cc1)OP(=O)(O)O |$_R1;;;;;_R2;;;;;;;;;;;;$| where * are connection points
 *
 * or * is square brackets like [*]N[C@H](C(=O)[*])Cc1ccc(cc1)OP(=O)(O)O |$_R1;;;;;_R2;;;;;;;;;;;;$|
 *
 * simple smiles with R notations like CCC[R1] or CCC(=O)[R2]
 *
 * simple smiles with ambiguety defined as [*:1] like CCC[*:1] or CCC(=O)[*:2]
*/
export function smiles2Monomer(smiles: string, polymerType?: PolymerType): MonomerWithoutSymbol | null {
  try {
    const isCxSmiles = cx_smiles_regexp.test(smiles);
    if (isCxSmiles) {
      // CXSMILES parsing
      const parts = smiles.split('|$');
      const molPart = parts[0].trim();
      const rGroupPart = parts[1];
      // make sure all R groups are captured
      const rGroupMatches = Array.from(rGroupPart.matchAll(/R(\d+)/g));
      const starsInMolecule = Array.from(molPart.matchAll(/(\*)/g));
      if (rGroupMatches.length !== starsInMolecule.length)
        return null; // make sure that number of R groups and stars are the same
      // remove brackets from stars if any
      let cleanMol = molPart.replaceAll(/\[\*\]/g, '*');
      // speaking in terms of consecutiveness, R groups in definition and stars in Smiles will be in the same order
      // so we can just iterate through them
      const rGroupNumbers = rGroupMatches.map((m) => m[1]); // numbers as strings
      for (let i = 0; i < rGroupNumbers.length; i++) {
        const rNum = rGroupNumbers[i];
        // replace first matched star with R group
        cleanMol = cleanMol.replace('*', `[R${rNum}]`);
      }
      return getMonomerFromRSmiles(cleanMol, polymerType);
    }

    // simple smiles parsing
    // to simplify, replace all ambigous [*:1] with R1, etc
    let cleanSmiles = smiles;
    const ambigMatches = Array.from(smiles.matchAll(ambig_regexp));
    for (const match of ambigMatches) {
      const fullMatch = match[0];
      const rNum = match[1];
      cleanSmiles = cleanSmiles.replace(fullMatch, `[R${rNum}]`);
    }

    // make sure monomer has at least one R group
    if (rgroup_regexp.test(cleanSmiles))
      return getMonomerFromRSmiles(cleanSmiles, polymerType);
  } catch (e) {
    _package.logger.error(`smiles2Monomer: cannot parse SMILES: ${smiles}\n${e}`);
  }


  return null;
}
