import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {ConnectionList} from './connection-list';
import {HELM_ITEM_SEPARATOR, HELM_SECTION_SEPARATOR} from './const';
import {SimplePolymer} from './simple-polymer';
import {Bond} from './types';

export class Helm {
  constructor(private helmString: string) {
    const helmSections = this.helmString.split(HELM_SECTION_SEPARATOR);
    const simplePolymers = helmSections[0].split(HELM_ITEM_SEPARATOR);
    this.simplePolymers = simplePolymers
      .map((item) => new SimplePolymer(item));
    this.connectionList = new ConnectionList(helmSections[1]);
    this.bondData = this.getBondData();

    this.bondedRGroupsMap = this.getBondedRGroupsMap();
  }

  /** List of pairs for bonded monomers, monomers indexed globally (withing the
   * complex polymer scope) */
  readonly bondData: Bond[][];

  public readonly simplePolymers: SimplePolymer[];
  public readonly connectionList: ConnectionList;

  /** Maps global monomer index to r-group ids (starting from 1) participating
   * in connection */
  readonly bondedRGroupsMap: number[][];

  private getBondedRGroupsMap(): number[][] {
    const monomerCount = this.simplePolymers.map((sp) => sp.monomers.length)
      .reduce((a, b) => a + b, 0);
    const bondedRGroupsList: number[][] = Array.from({length: monomerCount}, () => []);
    this.bondData.forEach((bond) => {
      bond.forEach((bondPart) => {
        const monomerIdx = bondPart.monomerIdx;
        const rGroupId = bondPart.rGroupId;
        bondedRGroupsList[monomerIdx].push(rGroupId);
      });
    });

    return bondedRGroupsList;
  }

  toString() {
    return this.helmString;
  }

  getPolymerTypeByMonomerIdx(monomerGlobalIdx: number): HELM_POLYMER_TYPE {
    const simplePolymer = this.getSimplePolymerByMonomerIdx(monomerGlobalIdx);
    const polymerType = simplePolymer.polymerType;
    return polymerType as HELM_POLYMER_TYPE;
  }

  private getSimplePolymerByMonomerIdx(monomerGlobalIdx: number): SimplePolymer {
    const shifts = this.getMonomerIdxShifts();
    const shiftValues = Object.values(shifts);
    const lowerBound = shiftValues.sort((a, b) => a - b).find(
      (shift) => monomerGlobalIdx >= shift
    );
    if (lowerBound === undefined)
      throw new Error(`Cannot find simple polymer for monomer ${monomerGlobalIdx}`);
    const simplePolymerId = Object.keys(shifts).find((simplePolymerId) => shifts[simplePolymerId] === lowerBound)!;
    const simplePolymer = this.simplePolymers.find((simplePolymer) => simplePolymer.id === simplePolymerId)!;
    return simplePolymer;
  }

  private shiftBondMonomerIds(shift: number, bonds: Bond[][]): void {
    bonds.forEach((bond) => {
      bond.forEach((bondPart) => {
        bondPart.monomerIdx += shift;
      });
    });
  }

  private getMonomerIdxShifts(): { [simplePolymerId: string]: number } {
    const result: { [simplePolymerId: string]: number } = {};
    let shift = 0;
    this.simplePolymers.forEach((simplePolymer) => {
      result[simplePolymer.id] = shift;
      shift += simplePolymer.monomers.length;
    });
    return result;
  }

  private getBondData(): Bond[][] {
    const shifts = this.getMonomerIdxShifts();
    const result: Bond[][] = [];
    this.simplePolymers.forEach((simplePolymer) => {
      const bondData = simplePolymer.getBondData();
      const shift = shifts[simplePolymer.id];
      this.shiftBondMonomerIds(shift, bondData);
      result.push(...bondData);
    });
    const connectionData = this.connectionList.getConnectionData();
    connectionData.forEach((connection) => {
      const data: Bond[] = [];
      connection.forEach((connectionItem) => {
        const shift = shifts[connectionItem.polymerId];
        const bond = connectionItem.bond;
        bond.monomerIdx += shift;
        data.push(bond);
      });
      result.push(data);
    });
    return result;
  }
}

