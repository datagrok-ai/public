import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IPdbHelper} from '@datagrok-libraries/bio';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';

import * as NGL from 'NGL';

export class PdbHelper implements IPdbHelper {

  //private stage: NGL.Stage;

  constructor() {
    //this.stage = new NGL.Stage();
  }

  async pdbToDf(pdbStr: string, name: string): Promise<DG.DataFrame> {
    const pdbBlob: Blob = new Blob([pdbStr], {type: 'text/plain'});
    //await this.stage.loadFile(pdbBlob, {ext: 'pdb'});

    //https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    const resDf: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('ResSeqNum', new Int32Array(0)), // Residue sequence number 23-26
    ]);
    const pdbObj: object = {type: 'pdb', message: 'Not implemented'};
    resDf.setTag(pdbTAGS.PDB, pdbStr);
    resDf.temp.set(pdbTAGS.PDB, pdbObj);

    return resDf;
  }
}