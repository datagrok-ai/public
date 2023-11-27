import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {delay} from '@datagrok-libraries/utils/src/test';
import {Molecule3DUnits} from '@datagrok-libraries/bio/src/molecule-3d/molecule-3d-units-handler';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb/index';

import {Pdbqt} from './pdbqt-parser';
import {IMPORT} from '../../consts-import';
import {buildDataJson, PROPS as bsPROPS} from '../../viewers/molstar-viewer/molstar-viewer';
import {Base64} from 'js-base64';
import {viewBiostructure} from '../../viewers/view-preview';


export async function importPdbqtUI(fileContent: string, test: boolean = false): Promise<DG.DataFrame[]> {
  const data: Pdbqt = Pdbqt.parse(fileContent);

  if (data.models.length > 0) {
    const molColName = IMPORT[Molecule3DUnits.pdbqt].molColName;
    const df: DG.DataFrame = data.toDataFrame(molColName);
    if (data.target) {
      const targetPdbStr = data.target.toPdb();
      df.setTag(pdbTAGS.PDB, targetPdbStr);
      return [df];
    } else {
      if (test) return [df];

      const dataFileProp = DG.Property.fromOptions({name: 'dataFile', caption: 'Data file', type: 'file'});
      const dataFileInput = DG.InputBase.forProperty(dataFileProp);

      const view = grok.shell.addTableView(df);
      const openTargetDlg = ui.dialog({title: 'Open file', showHeader: true})
        .add(ui.divText('Docking target structure required to display ligand poses from pdbqt data.'))
        .add(ui.divV([ui.inputs([dataFileInput])]))
        .onOK(async () => {
          await delay(100); /* to fill DG.FileInfo.data */
          const dataFi: DG.FileInfo = dataFileInput.value;
          const dataA: Uint8Array = dataFi.data ? dataFi.data /* User's file*/ : await dataFi.readAsBytes()/* Shares */;
          const viewer = await df.plot.fromType('Biostructure', {
            [bsPROPS.dataJson]: buildDataJson(dataA, dataFi.extension),
            [bsPROPS.ligandColumnName]: molColName,
          });
          view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Docking target biostructure', 0.4);
        })
        .show({centerAt: view.root});
    }
    return []; // return empty list, so open views required
  } else if (data.target) {
    const pdbStr = data.target!.toPdb();
    await viewBiostructure(pdbStr);
    return [];
  } else {
    grok.shell.warning('Content is not supported.');
    return [];
  }
}
