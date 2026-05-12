// AutoDock pose → receptor lookup. Mirrors the helpers in the Docking
// package's `utils/utils.ts` so BSV can resolve receptors locally for the
// AutoDock-pose context panel (`dockingInteractionsWidget` in `package.ts`)
// without having to call back into the Docking package.
//
// The receptor lookup convention is:
//   1. Parse the `REMARK 1 receptor.<name> J.` line from the pose to get
//      the receptor name (set by AutoDock when it writes a pose file).
//   2. Look in `System:AppData/Docking/targets/<receptorName>/` for a
//      `.pdbqt` or `.pdb` file (deposited by the Docking package's own
//      `files/targets/` content on install).
//   3. If not found locally (the user has BSV without Docking installed,
//      or the targets folder was wiped), fall back to fetching the PDB by
//      ID from RCSB.
//
// The `System:AppData/Docking/targets` path is hard-coded for now because
// that's where AutoDock poses reference receptors — making it configurable
// would only matter if a non-AutoDock docking tool emitted poses with the
// same `REMARK 1 receptor.` convention. Not a realistic scenario today.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {BiostructureData} from '@datagrok-libraries/bio/src/pdb/types';

const DOCKING_TARGET_PATH = 'System:AppData/Docking/targets';

/** AutoDock writes the receptor name into the pose file as a REMARK line:
 *  `REMARK  1  receptor. <name> J.`. The cheap regex match below is used
 *  as a per-row pre-screen in the batch loop (so we don't hit `getReceptorFile`
 *  with its 3 file-system RPCs on rows that can't possibly resolve). */
export const RECEPTOR_NAME_RE = /REMARK\s+1\s+receptor\.\s+(\S+)\s+J\./;

/** Resolve the receptor for an AutoDock pose. Looks up the receptor name
 *  from the pose's REMARK header, then either reads the matching file
 *  from `System:AppData/Docking/targets/<name>/` or falls back to fetching
 *  from RCSB by PDB ID. */
export async function getReceptorData(pdb: string): Promise<BiostructureData> {
  const match = pdb.match(RECEPTOR_NAME_RE);
  const receptorName = match ? match[1] : '';
  const receptor = await getReceptorFile(receptorName);
  return {
    binary: false,
    data: receptor
      ? await grok.dapi.files.readAsText(receptor)
      : await fetchPdbContent(receptorName.toUpperCase()),
    ext: receptor ? receptor.extension : 'pdb',
    options: {name: receptorName},
  } as BiostructureData;
}

async function getReceptorFile(receptorName: string): Promise<DG.FileInfo | undefined> {
  const receptorPath = `${DOCKING_TARGET_PATH}/${receptorName}`;
  const folderExists = await grok.dapi.files.exists(receptorPath);
  if (!folderExists) return undefined;
  const files = await grok.dapi.files.list(receptorPath);
  return files.find((file) => ['pdbqt', 'pdb'].includes(file.extension));
}

async function fetchPdbContent(pdbId: string, format: string = 'pdb'): Promise<string> {
  if (!pdbId) return '';
  try {
    const resp = await grok.dapi.fetchProxy(`https://files.rcsb.org/download/${pdbId}.${format}`);
    if (resp.ok) return await resp.text();
  } catch (_e) { /* RCSB unreachable — return empty so the caller can show an error */ }
  return '';
}
