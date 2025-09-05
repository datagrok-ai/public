import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { RevvityUser } from "./revvity-api";
import { openRevvityNode } from './view-utils';


export type RevvityLibrary = {
  id: string;
  name: string;
  types: RevvityType[];
}

export type RevvityType = {
  name: string;
  count: number;
}

export let libraries: RevvityLibrary[] | undefined = undefined;

export async function getRevvityLibraries(): Promise<RevvityLibrary[]> {
    if (!libraries) {
        const librariesStr = await grok.functions.call('RevvitySignalsLink:getLibraries');
        libraries = JSON.parse(librariesStr);
    }
    return libraries!;
}

export async function createInitialSatistics(statsDiv: HTMLDivElement) {

  ui.setUpdateIndicator(statsDiv, true, 'Loading statistics...');
  const libs = await getRevvityLibraries();
  const libObjForTable: any[] = [];
  for (const lib of libs) {
    for (const libType of lib.types)
      libObjForTable.push({ libName: lib.name, libType: libType.name, count: libType.count });
  }


  const table = ui.table(libObjForTable, (lib) => ([
    lib.libName ?? '',
    lib.libType ?? '',
    ui.link(lib.count, () => {
      const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
      node.expanded = true;
      openRevvityNode(node, lib.libName, lib.libType);
    }),
  ]),
    ['Library', 'Type', 'Count']);

  statsDiv.append(table);
  ui.setUpdateIndicator(statsDiv, false);
}