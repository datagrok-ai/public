import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {RGroupHandler} from './r-group-handler';

export abstract class MolfileWrapper {
  protected monomerSymbol: string;

  protected atoms: MolfileAtoms;
  protected bonds: MolfileBonds;
  protected rGroups: RGroupHandler;

  abstract deleteBondLineWithSpecifiedRGroup(rGroupId: number): void;
  abstract shiftCoordinates(shift: {x: number, y: number}): void;
  abstract rotateCoordinates(angle: number): void;
  abstract getBondLines(): string[];
  abstract getAtomLines(): string[];
  abstract removeRGroups(rGroupIds: number[]): void;
  abstract replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void;
  abstract getAttachmentAtomByRGroupId(rgroupId: number): number;
  abstract shiftBonds(shift: number): void;
  abstract capRGroups(capGroupElements: string[]): void;

  protected shiftR1GroupToOrigin(): void {
    const r1Idx = this.rGroups.getAtomicIdx(1);
    if (r1Idx === null)
      throw new Error(`Cannot find R1 group for monomer ${this.monomerSymbol}`);
    const {x, y} = this.atoms.atomCoordinates[r1Idx];
    this.atoms.shift({x: -x, y: -y});
  }

  protected alignR2AlongX(): void {
    const r2Idx = this.rGroups.getAtomicIdx(2);
    if (r2Idx === null)
      throw new Error(`Cannot find R2 group for monomer ${this.monomerSymbol}`);
    const r2Coordinates = this.atoms.atomCoordinates[r2Idx];
    const tan = r2Coordinates.y / r2Coordinates.x;
    const angle = Math.atan(tan);
    if (isNaN(angle))
      throw new Error(`Cannot calculate angle for R2 group for monomer ${this.monomerSymbol}`);
    this.rotateCoordinates(-angle);
  }

  protected shiftMonomerToDefaultPosition(): void {
    this.shiftR1GroupToOrigin();
    const r2Idx = this.rGroups.getAtomicIdx(2);
    if (r2Idx !== null)
      this.alignR2AlongX();
  }
}

