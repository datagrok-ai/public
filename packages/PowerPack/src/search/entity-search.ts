/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from '../widgets/web-widget';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {list} from 'datagrok-api/ui';

export async function functionSearch(s: string): Promise<any[]> {
  s = s.toLowerCase();
  return DG.Func.find()
    .filter((value) =>
      value.friendlyName.toLowerCase().includes(s) ||
      value.name.toLowerCase().includes(s) ||
      value.description?.toLowerCase()?.includes(s));
}

export async function appSearch(s: string): Promise<any[]> {
  s = s.toLowerCase();
  return DG.Func.find({tags: ['app']})
    .filter((value) =>
      value.friendlyName.toLowerCase().includes(s) ||
      value.name.toLowerCase().includes(s) ||
          value.description?.toLowerCase()?.includes(s));
}

export async function scriptsSearch(s: string): Promise<any[]> {
  s = s.toLowerCase();
  return (await grok.dapi.scripts.filter(s).list());
}

export async function usersSearch(s: string): Promise<any[]> {
  s = s.toLowerCase();
  return (await grok.dapi.users.filter(s).list());
}

function iframe(src: string): DG.Widget {
  return new WebWidget({
    src: src,
    width: '100%',
    height: '500px'
  });
}

export async function pubChemSearch(s: string): Promise<any | null> {
  return s !== 'aspirin' ?
    null :
    iframe( 'https://pubchem.ncbi.nlm.nih.gov/compound/aspirin#section=3D-Conformer&embed=true');
}

export async function pdbSearch(s: string): Promise<HTMLElement|null> {
  let prom :HTMLElement|null=null;
  if (s.length!=4) return prom;
  const up=s.toUpperCase();
  if (!/^[A-Z0-9]+$/.test(up)) return prom;

  // return (s.length == 4) ? iframe(`https://www.rcsb.org/3d-view/${s}`) : null;
  const response=await grok.dapi.fetchProxy('https://data.rcsb.org/rest/v1/holdings/status/'+up);
  if (response.ok) {
    const json=await response.json();
    if (json['rcsb_repository_holdings_combined_entry_container_identifiers'].update_id) {
      console.log(up + ' matches pdb!');
      prom = ui.link('PDB: '+up, 'https://www.rcsb.org/structure/' + up);
    } else {
      return prom;
    }
  }

  return prom;
}

export async function wikiSearch(s: string): Promise<DG.Widget | null> {
  return (s.toLowerCase().startsWith('wiki:')) ?
    iframe(`https://en.m.wikipedia.org/wiki/${s.substring(5).trim()}`) :
    null;
}
