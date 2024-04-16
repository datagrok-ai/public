import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {V2K_CONST, V3K_CONST} from './const';
import {Helm} from './helm';
import {MonomerWrapper} from './monomer-wrapper';

export class Polymer {
  constructor(helmString: string, private rdKitModule: RDModule) {
    this.helm = new Helm(helmString);
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;

  addMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const monomerWrapper = new MonomerWrapper(monomerSymbol, monomerIdx, this.helm, shift, this.rdKitModule);

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

  compileToMolfile(): string {
    const molfileHeader = '\nDatagrok\n';
    const atomLines: string[] = [];
    const bondLines: string[] = [];

    const atomNumberShifts = this.getAtomNumberShifts();
    this.monomerWrappers.forEach((monomerWrapper, idx) => {
      monomerWrapper.shiftBonds(atomNumberShifts[idx]);
    });

    this.restoreBondsBetweenMonomers();

    this.monomerWrappers.forEach((monomerWrapper) => {
      atomLines.push(...monomerWrapper.getAtomLines());
      bondLines.push(...monomerWrapper.getBondLines());
    });

    const atomCount = atomLines.length;
    const bondCount = bondLines.length;

    const header = this.getV3KHeader(atomCount, bondCount);
    const atomBlock = this.getV3KAtomBlock(atomLines);
    const bondBlock = this.getV3KBondBlock(bondLines);
    const molfileEnd = V3K_CONST.END_CTAB + '\n' + V3K_CONST.END;
    const blockList = [header, atomBlock, bondBlock, molfileEnd];
    const molfile = blockList.join('\n');
    return molfile;
  }

  private getV3KHeader(atomCount: number, bondCount: number): string {
    const countsLine = `${V3K_CONST.COUNTS_LINE_START}${atomCount} ${bondCount}${V3K_CONST.COUNTS_LINE_END}`;
    return `${V3K_CONST.HEADER}\n${V3K_CONST.BEGIN_CTAB}\n${countsLine}`;
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
