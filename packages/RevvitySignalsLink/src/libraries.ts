import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { RevvityUser } from "./revvity-api";


export type RevvityLibrary = {
  id: string;
  name: string;
  types: string[];
}

export let libraries: RevvityLibrary[] | undefined = undefined;

export async function getRevvityLibraries(): Promise<RevvityLibrary[]> {
    if (!libraries) {
        const librariesStr = await grok.functions.call('RevvitySignalsLink:getLibraries');
        libraries = JSON.parse(librariesStr);
    }
    return libraries!;
}