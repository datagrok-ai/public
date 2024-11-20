import {IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {Linkage} from './pt-misc';
import {Rules} from './pt-rules';
import {HelmMol} from '@datagrok-libraries/bio/src/helm/types';
import {fromObjectsToHelm, handleDuplicated, handleLinkRules,
  handleReactionRules, parseHelm, parseSeparator} from './pt-tools-parse';
import {getHelmMol, helmMolToNotation} from './pt-tools-helmmol';


export class Chain {
  linkages: Linkage[];
  monomers: string[][];
  mol: HelmMol;

  underRules: boolean = false;
  linkagesUnderRules: Linkage[];
  monomersUnderRules: string[][];
  molUnderRules: HelmMol;

  constructor(monomers: string[][], linkages: Linkage[], protected helmHelper: IHelmHelper) {
    this.linkages = linkages;
    this.monomers = monomers;
    this.mol = getHelmMol(linkages, monomers, helmHelper);
  }

  /** Parse from separator or sequence notation (template)  */
  static fromSeparator(sequence: string, helmHelper: IHelmHelper): Chain {
    const [linkages, mainFragments] = parseSeparator(sequence);
    const chain = new Chain(mainFragments, linkages, helmHelper);
    return chain;
  }

  /** Parse harmonized sequence (template) from pseudo helm */
  static fromHelm(sequence: string, helmHelper: IHelmHelper) : Chain {
    const [linkages, mainFragments] = parseHelm(sequence);
    const chain = new Chain(mainFragments, linkages, helmHelper);
    return chain;
  }

  getNotation(): string {
    return helmMolToNotation(this.mol);
  }

  getHelm(): string {
    if (this.underRules)
      return fromObjectsToHelm(this.linkagesUnderRules, this.monomersUnderRules);
    else
      return fromObjectsToHelm(this.linkages, this.monomers);
  }

  /** Get macromolecule from harmonized sequence (template) */
  applyRules(rules: Rules): void {
    const sequence = this.getNotation();

    const [linkages, mainFragments] = handleDuplicated(sequence, rules);
    const monomers = new Array<Array<string>>(mainFragments.length);
    handleLinkRules(mainFragments, monomers, linkages, rules);
    handleReactionRules(monomers, linkages, rules);

    this.underRules = true;
    this.linkagesUnderRules = linkages;
    this.monomersUnderRules = monomers;
    this.molUnderRules = getHelmMol(linkages, monomers, this.helmHelper);
  }

  public check(throwError: boolean = false): string[] {
    const errors: string[] = [];

    const chainsMonomerCount = this.monomers.map((ch) => ch.length).reduce((acc, curr) => acc + curr, 0);
    if (this.mol.atoms.length !== chainsMonomerCount) {
      errors.push(`The mol atoms count ${this.mol.atoms.length} does not match ` +
        `the total number ${chainsMonomerCount} of chains' monomers.`);
    }

    const internalBondsCount = this.monomers.map((ch) => ch.length - 1).reduce((acc, curr) => acc + curr, 0);
    const chainsBondCount = internalBondsCount + this.linkages.length;
    if (this.mol.bonds.length !== chainsBondCount) {
      errors.push(`The mol bonds count ${this.mol.bonds.length} does not match ` +
        `the total number ${chainsBondCount} in- and inter-chain linkages.`);
    }

    let counter: number = 0;
    for (let spIdx = 0; spIdx < this.monomers.length; ++spIdx) {
      const chain = this.monomers[spIdx];
      for (let mIntIdx = 0; mIntIdx < chain.length; ++mIntIdx) {
        try {
          const m = chain[mIntIdx];
          const a = this.mol.atoms[counter];
          if (a.bio!.continuousId !== counter)
            errors.push(`Atom #${counter} has incorrect .bio.continuousId: ${a.bio!.continuousId}.`);
          if (a.elem !== m)
            errors.push(`Atom #${counter} elem: '${a.elem}' does not match chain monomer: '${m}'.`);
        } finally { counter++; }
      }
    }
    if (throwError && errors.length > 0)
      throw new Error(`Chain errors:\n${errors.map((e) => `  ${e}`).join('\n')}`);
    return errors;
  }
}
