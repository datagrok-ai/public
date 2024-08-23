import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {RGroupHandler} from './r-group-handler';

export abstract class MolfileWrapper {
  constructor(protected monomerSymbol: string) { }

  protected atoms: MolfileAtoms;
  protected bonds: MolfileBonds;
  protected rGroups: RGroupHandler;

  public get atomCount(): number { return this.atoms.count; }

  public get bondCount(): number { return this.bonds.count; }

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

  deleteBondLineWithSpecifiedRGroup(rGroupId: number): void {
    this.rGroups.deleteBondLineWithSpecifiedRGroup(rGroupId);
  }

  shiftCoordinates(shift: { x: number, y: number }): void {
    this.atoms.shift(shift);
  }

  rotateCoordinates(angle: number): void {
    this.atoms.rotate(angle);
  }

  getBondLines(): string[] {
    return this.bonds.getBondLines();
  }

  getAtomLines(): string[] {
    return this.atoms.atomLines;
  }

  removeRGroups(rGroupIds: number[]): void {
    this.rGroups.removeRGroups(rGroupIds);
  }

  replaceRGroupWithAttachmentAtom(rGroupId: number, externalAtom: number): void {
    this.rGroups.replaceRGroupWithAttachmentAtom(rGroupId, externalAtom);
  }

  getAttachmentAtomByRGroupId(rgroupId: number): number {
    return this.rGroups.getAttachmentAtomIdByRGroupId(rgroupId);
  }

  shiftBonds(shift: number): void {
    this.bonds.shift(shift);
  }

  capRGroups(capGroupElements: string[]): void {
    this.rGroups.capRGroups(capGroupElements);
  }
}

