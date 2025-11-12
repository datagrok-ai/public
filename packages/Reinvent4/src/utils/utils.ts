import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import JSZip from 'jszip';

export async function zipFolder(folder: DG.FileInfo[]): Promise<Blob> {
  const zip = new JSZip();
  folder.forEach(file => {
    zip.file(file.name, file.readAsBytes());
  });

  const zipBlob = await zip.generateAsync({ type: 'blob' });
  return zipBlob;
}