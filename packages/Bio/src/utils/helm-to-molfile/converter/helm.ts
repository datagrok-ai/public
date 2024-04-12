import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {ConnectionList} from './connection-list';
import {HELM_ITEM_SEPARATOR, HELM_SECTION_SEPARATOR} from './const';
import {SimplePolymer} from './simple-polymer';
import {Bond} from './types';

export class Helm {
  constructor(private helm: string) {
    const helmSections = this.helm.split(HELM_SECTION_SEPARATOR);
    const simplePolymers = helmSections[0].split(HELM_ITEM_SEPARATOR);
    this.simplePolymers = simplePolymers
      .map((item) => new SimplePolymer(item));
    if (helmSections[1] !== '')
      this.connectionList = new ConnectionList(helmSections[1]);
    this.bondData = this.getBondData();
  }

  /** List of pairs for bonded monomers, monomers indexed globally (withing the
   * complex polymer scope) */
  readonly bondData: Bond[][];

  private simplePolymers: SimplePolymer[];
  private connectionList?: ConnectionList;

  toString() {
    return this.helm;
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

  private getMonomerIdxShifts(): {[simplePolymerId: string]: number} {
    const result: {[simplePolymerId: string]: number} = {};
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
    if (this.connectionList) {
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
    }
    return result;
  }
}

