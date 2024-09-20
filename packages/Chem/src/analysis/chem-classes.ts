import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getMmpFrags} from './molecular-matched-pairs/mmp-analysis/mmpa-fragments';

// export interface IMmpFragmentsResult {
//   frags: [string, string][][];
//   smiles: string[];
// }

const CLASS_CUTOFF = 0.8;

export async function getChemClasses(molecules: DG.Column): Promise<string[]> {
  const frags = await getMmpFrags(molecules.toList());

  let fragmentCounter = 0;
  //const encoded: [number, number][][] = new Array(frags.frags.length)
  //.fill(null).map((_, i) => new Array(frags.frags[i].length).fill(null).map((_) => [0, 0]));
  const encoded: number[][] = new Array(frags.frags.length)
    .fill(null).map((_, i) => new Array(frags.frags[i].length).fill(null).map((_) => 0));

  const fragmentMap: Record<string, number> = {};
  const fragmentSumsMap: Record<string, number> = {};
  for (let i = 0; i < frags.frags.length; i++) {
    for (let j = 0; j < frags.frags[i].length; j++) {
      if (!fragmentMap[frags.frags[i][j][0]]) {
        fragmentMap[frags.frags[i][j][0]] = fragmentCounter;
        fragmentCounter++;
        fragmentSumsMap[frags.frags[i][j][0]] = 1;
      } else
        fragmentSumsMap[frags.frags[i][j][0]]++;
      // if (!fragmentMap[frags.frags[i][j][1]]) {
      //   fragmentMap[frags.frags[i][j][1]] = fragmentCounter;
      //   fragmentCounter++;
      // }

      encoded[i][j]/*[0]*/= fragmentMap[frags.frags[i][j][0]];
      //[i][j][1] = fragmentMap[frags.frags[i][j][1]];
    }
  }

  // const maxFragmentIndex = fragmentCounter;
  // const fragIdToFragName = new Array<string>(maxFragmentIndex);
  // Object.entries(fragmentMap).forEach(([key, val]) => {
  //   fragIdToFragName[val] = key;
  // });

  const fragsUnique = Object.keys(fragmentSumsMap).sort((a: string, b: string) => a.length > b.length ? -1 : 1);
  let frag = '';
  const cuttoff = CLASS_CUTOFF * molecules.length;

  for (let i = 0; i < fragsUnique.length; i++) {
    if (fragmentSumsMap[fragsUnique[i]] >= cuttoff) {
      frag = fragsUnique[i];
      break;
    }
  }

  const res = new Array<string>(molecules.length).fill(frag);

  return res;
}
