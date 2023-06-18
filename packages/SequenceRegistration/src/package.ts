/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {engageViewForOligoSdFileUI} from './utils/registration';
import {DBLoaderBase, DataLoaderDB} from './utils/database-loader';

class Package extends DG.Package {
  private _dbLoader?: DBLoaderBase;

  get dbLoader(): DBLoaderBase {
    if (!this._dbLoader)
      throw new Error('ST: dataLoader is not initialized');
    return this._dbLoader!;
  };

  public async initDBLoader(): Promise<void> {
    if (this._dbLoader === undefined) {
      const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
        'Initializing Sequence Registration data loader ...');
      try {
        const dl = new DataLoaderDB();
        await dl.init();
        this._dbLoader = dl;
      } catch (err: any) {
        const errMsg = err.hasOwnProperty('message') ? err.message : err.toString();
        grok.shell.error(errMsg);
        throw new Error('Initializing Sequence Registration data loader error: ' + errMsg);
      } finally {
        pi.close();
      }
    }
  }
}

export const _package = new Package();

//tags: init
export async function initSequenceRegistration(): Promise<void> {
  _package.initDBLoader().then(() => {}); // Do not wait for lists loaded from the database
}

//name: engageViewForOligoSdFile()
//description: Function to modify the view for oligo batch registration
export async function engageViewForOligoSdFile(view: DG.TableView): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
    'Processing table for Sequence Registration...');
  await engageViewForOligoSdFileUI(view);
  pi.close();
}
