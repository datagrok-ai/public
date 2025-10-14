import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {HELM_ITEM_SEPARATOR} from './const';
import {Bond} from './types';

export class ConnectionList {
  constructor(connectionList: string) {
    const splitted = connectionList.split(HELM_ITEM_SEPARATOR).filter((ci) => ci);
    splitted.forEach((connectionItem: string) => this.validateConnectionItem(connectionItem));
    this.connectionItems = splitted;
  }

  public connectionItems: string[];

  private validateConnectionItem(connectionItem: string): void {
    const allowedType = `(${HELM_POLYMER_TYPE.PEPTIDE}|${HELM_POLYMER_TYPE.RNA})`;
    const regex = new RegExp(`${allowedType}[0-9]+,${allowedType}[0-9]+,[0-9]+:R[0-9]+-[0-9]+:R[0-9]+`, 'g');
    if (!connectionItem.match(regex))
      throw new Error(`Cannot parse connection item from ${connectionItem}`);
  }

  getConnectionData(): { polymerId: string, bond: Bond }[][] {
    const result: { polymerId: string, bond: Bond }[][] = [];
    this.connectionItems.forEach((connectionItem: string) => {
      const pair: { polymerId: string, bond: Bond }[] = [];
      const splitted = connectionItem.split(',');
      splitted[2].split('-').forEach((item, idx) => {
        const polymerId = splitted[idx];
        const data = item.split(':');
        // WARNING: monomer idx starts from 0
        const monomerIdx = parseInt(data[0]) - 1;
        const rGroupId = parseInt(data[1].slice(1));
        const bondData = {monomerIdx, rGroupId};
        pair.push({polymerId, bond: bondData});
      });
      result.push(pair);
    });
    return result;
  }
}

