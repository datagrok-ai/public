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

export async function createInitialSatistics(statsDiv: HTMLElement, libName?: string) {

  ui.setUpdateIndicator(statsDiv, true, 'Loading statistics...');
  const libs = await getRevvityLibraries();
  const libObjForTable: any[] = [];
  for (const lib of libs) {
    if (libName && lib.name !== libName)
      continue;
    for (const libType of lib.types)
      libObjForTable.push({ libName: lib.name, libType: libType.name, count: libType.count });
  }

  const statsElement = createLibsStatsTable(libObjForTable, libName);
  statsDiv.append(statsElement);
  ui.setUpdateIndicator(statsDiv, false);
}

function createLibsStatsTable(libObjForTable: any[], libName?: string): HTMLElement {
  const output = libName ? ['Type', 'Count'] : ['Library', 'Type', 'Count'];

  const table = ui.table(libObjForTable, (lib) => {
    const arr = [
      lib.libName ?? '',
      lib.libType ?? '',
      ui.link(lib.count, () => {
        const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
        node.expanded = true;
        openRevvityNode(node, [lib.libName], lib.libType, lib.libName, lib.libType);
      }),
    ];
    if (libName)
      arr.splice(0, 1);

    return arr;
  },
    output);

  return libName ? ui.div([ui.h3(libName), table]) : table;
}