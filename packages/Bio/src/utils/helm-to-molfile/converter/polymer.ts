import wu from 'wu';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {V3K_CONST} from '@datagrok-libraries/chem-meta/src/formats/molfile-const';
import {HelmTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {IMonomerLib, IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GapOriginals} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {MolfileWithMap, MonomerMap} from '@datagrok-libraries/bio/src/monomer-works/types';

import {Helm} from './helm';
import {MonomerWrapper} from './monomer-wrapper';

export class Polymer {
  constructor(
    helmString: string,
    private readonly rdKitModule: RDModule,
    private readonly monomerLib: IMonomerLibBase
  ) {
    this.helm = new Helm(helmString);
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;

  addMonomer(
    monomerSymbol: string,
    helmMonomerIdx: number,
    shift: { x: number, y: number },
  ): void {
    if (monomerSymbol === GapOriginals[NOTATION.HELM]) return;

    const molMonomerIdx: number = this.monomerWrappers.length;
    const monomerWrapper = new MonomerWrapper(
      monomerSymbol, molMonomerIdx, this.helm, shift, this.rdKitModule, this.monomerLib);

    this.monomerWrappers.push(monomerWrapper);
  }

  private getAtomNumberShifts(): number[] {
    const atomNumberShifts: number[] = [];
    let shift = 0;
    this.monomerWrappers.forEach((monomerWrapper) => {
      atomNumberShifts.push(shift);
      shift += monomerWrapper.getAtomLines().length;
    });
    return atomNumberShifts;
  }

  private restoreBondsBetweenMonomers(): void {
    this.helm.bondData.forEach((bond) => {
      const monomerIdx = bond.map((bondPart) => bondPart.monomerIdx);
      const rGroupId = bond.map((bondPart) => bondPart.rGroupId);
      const monomer = monomerIdx.map((idx) => this.monomerWrappers[idx]);

      const attachmentAtom = monomer[1].getAttachmentAtomByRGroupId(rGroupId[1]);
      monomer[0].replaceRGroupWithAttachmentAtom(rGroupId[0], attachmentAtom);
      monomer[1].deleteBondLineWithSpecifiedRGroup(rGroupId[1]);
    });
  }

  compileToMolfile(): MolfileWithMap {
    const atomLines: string[] = [];
    const bondLines: string[] = [];

    const atomNumberShifts = this.getAtomNumberShifts();
    this.monomerWrappers.forEach((monomerWrapper, idx) => {
      monomerWrapper.shiftBonds(atomNumberShifts[idx]);
    });

    this.restoreBondsBetweenMonomers();

    const monomers: MonomerMap = new MonomerMap();
    for (const [mw, mwI] of wu.enumerate(this.monomerWrappers)) {
      const mAtomFirst = atomLines.length;
      const mBondFirst = bondLines.length;

      atomLines.push(...mw.getAtomLines());
      bondLines.push(...mw.getBondLines());

      const polymerType = this.helm.getPolymerTypeByMonomerIdx(mwI);
      const biotype = polymerType == PolymerTypes.RNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
      monomers.set(mwI, {
        biotype: biotype,
        symbol: mw.monomerSymbol,
        atoms: wu.count(mAtomFirst).take(mw.atomCount).toArray(),
        bonds: wu.count(mBondFirst).take(mw.bondCount).toArray(),
      });
    }

    const atomCount = atomLines.length;
    const bondCount = bondLines.length;

    const header = this.getV3KHeader(atomCount, bondCount);
    const atomBlock = this.getV3KAtomBlock(atomLines);
    const bondBlock = this.getV3KBondBlock(bondLines);
    const molfileEnd = V3K_CONST.END_CTAB + '\n' + V3K_CONST.END;
    const blockList = [header, atomBlock, bondBlock, molfileEnd];
    const molfile = blockList.join('\n');
    return {molfile, monomers};
  }

  private getV3KHeader(atomCount: number, bondCount: number): string {
    const countsLine = `${V3K_CONST.COUNTS_LINE_START}${atomCount} ${bondCount} 0 0 1`;
    return `${V3K_CONST.DUMMY_HEADER}\n${V3K_CONST.BEGIN_CTAB}\n${countsLine}`;
  }

  private getV3KAtomBlock(atomLines: string[]): string {
    const regex = /^(M  V30 )(\d+)( .*)$/;
    const newAtomLines = atomLines.map((line, idx) => {
      const atomIndex = idx + 1;
      return line.replace(regex, (match, p1, p2, p3) => { return p1 + atomIndex + p3; });
    });

    const atomBlock = [V3K_CONST.BEGIN_ATOM_BLOCK, ...newAtomLines, V3K_CONST.END_ATOM_BLOCK];
    return atomBlock.join('\n');
  }

  private getV3KBondBlock(bondLines: string[]): string {
    const regex = /^(M  V30 )(\d+)( .*)$/;
    const newBondLines = bondLines.map((line, idx) => {
      const atomIndex = idx + 1;
      return line.replace(regex, (match, p1, p2, p3) => { return p1 + atomIndex + p3; });
    });
    const bondBlock = [V3K_CONST.BEGIN_BOND_BLOCK, ...newBondLines, V3K_CONST.END_BOND_BLOCK];
    return bondBlock.join('\n');
  }
}
