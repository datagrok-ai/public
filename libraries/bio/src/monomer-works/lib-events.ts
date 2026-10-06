/* The platform-wide word that the monomer libraries have (re)loaded — fired by Bio's library
   manager after every load, for any package (or a test) that renders or validates monomers and has
   no library object of its own to subscribe to. Carried on `grok.events` custom events, so a
   listener needs neither Bio nor this library at hand: `grok.events.onCustomEvent('bio-monomer-lib-loaded')`. */
import * as grok from 'datagrok-api/grok';
import {Observable} from 'rxjs';

export const MONOMER_LIB_LOADED_EVENT = 'bio-monomer-lib-loaded';

export interface MonomerLibLoadedArgs {
  /** The sources loaded into the monomer library, the file names for the files provider. */
  libraries: string[];
  /** Whether the library was cleared before the load (a settings change) or merged into. */
  reload: boolean;
}

export function onMonomerLibLoaded(): Observable<MonomerLibLoadedArgs> {
  return grok.events.onCustomEvent(MONOMER_LIB_LOADED_EVENT);
}

export function fireMonomerLibLoaded(args: MonomerLibLoadedArgs): void {
  grok.events.fireCustomEvent(MONOMER_LIB_LOADED_EVENT, args);
}
