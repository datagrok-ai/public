import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {MolfileAtoms} from './mol-atoms';
import {MolfileBonds} from './mol-bonds';
import {MolfileWrapper} from './mol-wrapper';
import {RGroupHandler} from './r-group-handler';

export class MolfileV2KWrapper extends MolfileWrapper {
  constructor(molfileV2K: string, protected monomerSymbol: string) {
    super();
    const molfileHandler = MolfileHandler.getInstance(molfileV2K);

    this.atoms = new MolfileAtoms(molfileHandler);
    this.bonds = new MolfileBonds(molfileHandler);

    this.rGroups = new RGroupHandler(molfileHandler, this.atoms, this.bonds);

    this.shiftMonomerToDefaultPosition();
  }

  protected atoms: MolfileAtoms;
  protected bonds: MolfileBonds;
  protected rGroups: RGroupHandler;

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

