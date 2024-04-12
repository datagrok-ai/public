import {V2K_CONST} from './const';
import {Helm} from './helm';
import {MonomerWrapper} from './monomer-wrapper';

export class Polymer {
  constructor(helm: string) {
    this.helm = new Helm(helm);

    this.bondedRGroupsMap = new Map<number, number[]>();
    this.helm.bondData.forEach((bond) => {
      bond.forEach((bondPart) => {
        const monomerIdx = bondPart.monomerIdx;
        const rGroupId = bondPart.rGroupId;
        if (!this.bondedRGroupsMap.get(monomerIdx))
          this.bondedRGroupsMap.set(monomerIdx, []);
        this.bondedRGroupsMap.get(monomerIdx)!.push(rGroupId);
      });
    });
  }

  private monomerWrappers: MonomerWrapper[] = [];
  private helm: Helm;
  /** Maps global monomer index to r-group ids (starting from 1) participating
   * in connection */
  private bondedRGroupsMap: Map<number, number[]>;

  addMonomer(
    monomerSymbol: string,
    monomerIdx: number,
    shift: {x: number, y: number},
  ): void {
    const polymerType = this.helm.getPolymerTypeByMonomerIdx(monomerIdx);
    const monomerWrapper = new MonomerWrapper(monomerSymbol, polymerType);
    monomerWrapper.shiftCoordinates(shift);

    this.monomerWrappers.push(monomerWrapper);
  }

  private removeRGroups(): void {
    this.monomerWrappers.forEach((monomerWrapper, monomerIdx) => {
      if (this.bondedRGroupsMap.has(monomerIdx))
        monomerWrapper.removeBondedRGroups(this.bondedRGroupsMap.get(monomerIdx)!);
      monomerWrapper.capTrailingRGroups();
    });
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

    this.removeRGroups();

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
