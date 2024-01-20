/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from '../widgets/web-widget';

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

export async function pubChemSearch(s: string): Promise<DG.Widget | null> {
  return s !== 'aspirin' ?
    null :
    iframe( 'https://pubchem.ncbi.nlm.nih.gov/compound/aspirin#section=3D-Conformer&embed=true');
}

export async function pdbSearch(s: string): Promise<DG.Widget | null> {
  return (s.length == 4) ? iframe(`https://www.rcsb.org/3d-view/${s}`) : null;
}

export async function wikiSearch(s: string): Promise<DG.Widget | null> {
  return (s.toLowerCase().startsWith('wiki:')) ?
    iframe(`https://en.m.wikipedia.org/wiki/${s.substring(5).trim()}`) :
    null;
}
