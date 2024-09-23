import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {deserialize, serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

declare global {
  var SNAPSHOTS_UPDATE_MODE: boolean;
}

const snapshotsPath = 'System:AppData/LibTests/snapshots/';

export async function snapshotCompare(actual: any, snapshotName: string) {
  // if (globalThis.SNAPSHOTS_UPDATE_MODE) {
  //   const data = serialize(actual, {useJsonDF: true, space: 2});
  //   const name = snapshotName + '.json';
  //   const blob = new Blob([data]);
  //   DG.Utils.download(name, blob);
  // } else {
  const data = await grok.dapi.files.readAsText(snapshotsPath + snapshotName + '.json');
  const expected = deserialize(data);
  expectDeepEqual(actual, expected);
  // }
}
