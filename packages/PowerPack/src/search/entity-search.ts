/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {WebWidget} from '../widgets/web-widget';
import {HelpObject} from './help-entity';
import {helpInfo} from './helpIndex.g';


export async function functionSearch(s: string): Promise<DG.Func[]> {
  s = s.toLowerCase().trim();
  return DG.Func.find()
    .filter((value) => (value.name.toLowerCase().includes(s) || value.description?.toLowerCase()?.includes(s)) &&
    (!(value instanceof DG.Script) && !(value instanceof DG.DataQuery)) && !value.hasTag('app'));
}

export async function appSearch(s: string): Promise<DG.Func[]> {
  s = s.toLowerCase().trim();
  return DG.Func.find({tags: ['app']}).filter((val) => val.name?.toLowerCase().includes(s) ||
    val.description?.toLowerCase().includes(s) || val.friendlyName?.toLowerCase().includes(s));
}

export async function scriptsSearch(s: string): Promise<DG.Func[]> {
  s = s.toLowerCase().trim();
  return DG.Func.find()
    .filter((value) => (value.name.toLowerCase().includes(s) || value.description?.toLowerCase()?.includes(s)) &&
    (value instanceof DG.Script));
}

export async function querySearch(s: string): Promise<DG.Func[]> {
  s = s.toLowerCase().trim();
  return DG.Func.find()
    .filter((value) => (value.name.toLowerCase().includes(s) || value.description?.toLowerCase()?.includes(s)) &&
    (value instanceof DG.DataQuery));
}

export async function connectionsSearch(s: string): Promise<DG.DataConnection[]> {
  s = s.toLowerCase().trim();
  return (await grok.dapi.connections.filter(s).list());
}
export async function usersSearch(s: string): Promise<DG.User[]> {
  s = s.toLowerCase().trim();
  return (await grok.dapi.users.filter(s).list());
}

export async function helpSearch(s: string): Promise<HelpObject[]> {
  s = s.toLowerCase().trim();
  const getMatchScore = (h: HelpObject): number => {
    let score = 0;
    const addMatches = (ss: string, multiplier: number): void => {
      if (h.title.toLowerCase().includes(ss)) score += 10 * multiplier;
      h.keywords?.filter((k) => k.toLowerCase().includes(ss))?.forEach(() => score += 5 * multiplier);
      if (h.helpURL.toLowerCase().includes(ss)) score += 3 * multiplier;
    };

    addMatches(s, 1);
    const words = s.split(' ').map((w) => w.trim()).filter((w) => !!w);
    if (words.length > 1)
      words.forEach((w) => addMatches(w, 0.5));

    return score;
  };

  const matches = helpInfo.map((h) => ({o: h, score: getMatchScore(h)}))
    .filter((h) => h.score > 0)
    .sort((a, b) => b.score - a.score)
    .map((h) => HelpObject.fromHelpInfo(h.o));
  return matches;
}

export async function groupsSearch(s: string): Promise<DG.Group[]> {
  s = s.toLowerCase().trim();
  return (await grok.dapi.groups.filter(s).list()).filter((g) => g.personal === false);
}

export async function dockerSearch(s: string): Promise<DG.DockerContainer[]> {
  s = s.toLowerCase().trim();
  return await grok.dapi.docker.dockerContainers.filter(s).list();
}

function iframe(src: string): DG.Widget {
  return new WebWidget({
    src: src,
    width: '100%',
    height: '500px',
  });
}

export async function pubChemSearch(s: string): Promise<DG.Widget | null> {
  return s !== 'aspirin' ?
    null :
    iframe( 'https://pubchem.ncbi.nlm.nih.gov/compound/aspirin#section=3D-Conformer&embed=true');
}

export async function pdbSearch(s: string): Promise<DG.Widget | null> {
  if (s.length != 4 || !/^[A-Z0-9]+$/.test(s.toUpperCase()))
    return null;

  const response= await grok.dapi.fetchProxy('https://data.rcsb.org/rest/v1/holdings/status/' + s);
  if (!response.ok)
    return null;

  const json = await response.json();
  if (json['rcsb_repository_holdings_combined_entry_container_identifiers']?.['update_id']) {
    return DG.Widget.fromRoot(ui.divH([
      ui.link('PDB: ' + s, 'https://www.rcsb.org/structure/' + s),
      ui.image(`https://cdn.rcsb.org/images/structures/${s.toLowerCase()}_assembly-1.jpeg`, 200, 150),
    ]));
  }

  return null;
}

export async function wikiSearch(s: string): Promise<DG.Widget | null> {
  return (s.toLowerCase().startsWith('wiki:')) ?
    iframe(`https://en.m.wikipedia.org/wiki/${s.substring(5).trim()}`) :
    null;
}

// NOTE: this is an easter egg :D
export async function denialSearch(s: string, w?: DG.Widget): Promise<DG.Widget | null> {
// This function is a joke :D
// in the era of web services for everything, god gave us a service for generating deniel messages
  s = s.toLowerCase().trim();
  const nos = ['deny', 'say no', 'generate denial', 'how do i say no politely', 'generate no reason', 'generate no'];
  if (!nos.some((no) => s.includes(no)))
    return null;
  try {
    const response = await grok.dapi.fetchProxy('https://naas.isalman.dev/no', {method: 'GET'});
    if (!response.ok)
      return null;
    const json = await response.json();
    if (json?.reason) {
      const widget = w ?? DG.Widget.fromRoot(ui.div([], {style: {display: 'flex'}}));
      Array.from(widget.root.children).forEach((c) => c.remove());
      const refreshIcon = ui.icons.sync(() => denialSearch(s, widget));
      refreshIcon.style.marginLeft = '10px';
      const content =
        ui.divH([ui.h1(json.reason, {style: {textAlign: 'center', margin: 'auto'}}),
          refreshIcon], {style: {margin: 'auto', alignItems: 'center'}});
      widget.root.appendChild(content);
      return widget;
    }
    return null;
  } catch (e) {
    console.error('Error fetching denial message:', e);
    return null;
  }
}
