import {V2K_CONST} from './const';
import {Helm} from './helm';
import {MonomerWrapper} from './monomer-wrapper';

export class Polymer {
  constructor(helmString: string) {
    this.helm = new Helm(helmString);
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;

  addMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const monomerWrapper = new MonomerWrapper(monomerSymbol, monomerIdx, this.helm, shift);

    this.monomerWrappers.push(monomerWrapper);
  }

  private getAtomNumberShifts(): number[] {
    let shift = 0;
    const atomNumberShifts = this.monomerWrappers.map(
      (monomerWrapper) => shift += monomerWrapper.getAtomLines().length
    );
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
    if (atomCount > V2K_CONST.MAX_ATOM_COUNT) {
      throw new Error(
        `Atom count in polymer ${this.helm.toString()} is ${atomCount} and exceeds ${V2K_CONST.MAX_ATOM_COUNT}`
      );
    }

    const bondCount = bondLines.length;
    const countsLine = `${
      atomCount.toString().padStart(3, ' ')
    }${
      bondCount.toString().padStart(3, ' ')
    }  0  0  1  0              0 V2000`;
    const molfileEnd = 'M  END\n';
    const newLineChar = '\n';
    const blockList = [molfileHeader, countsLine, atomLines.join(newLineChar), bondLines.join(newLineChar), molfileEnd];
    const molfile = blockList.join(newLineChar);
    return molfile;
  }
}
