/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {KetcherSketcher} from './ketcher';

export let _package = new DG.Package();

//name: ketcherSketcher
//tags: moleculeSketcher
//output: widget sketcher
export function ketcherSketcher() {
  return new KetcherSketcher();
}
