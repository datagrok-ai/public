import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {Helm} from '../utils/helm-to-molfile/converter/helm';

category('helm', () => {
  const tests = {
    'single-linear': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T.[meI]}$$$$'},
      tgt: {
        simplePolymers: [6],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    },
    'single-cyclized-C-2-2': {
      src: {helm: 'PEPTIDE1{R.F.C.Y.G.H.[GGaz].C.T.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-8,R3$$$'},
      tgt: {
        simplePolymers: [10],
        bondedRGroups: [1, 2, 3, 2, 2, 2, 2, 3, 2, 1],
      }
    },
    'single-cyclized-C-1-1': {
      src: {helm: 'PEPTIDE1{F.C.Y.G.H.[GGaz].C.[meI]}$PEPTIDE1,PEPTIDE1,2:R3-7,R3$$$'},
      tgt: {
        simplePolymers: [8],
        bondedRGroups: [1, 3, 2, 2, 2, 1, 3, 1],
      }
    },
    'single-cyclized-C-0-0': {
      src: {helm: 'PEPTIDE1{C.Y.G.H.[GGaz].C}$PEPTIDE1,PEPTIDE1,2:R3-7,R3$$$'},
      tgt: {
        simplePolymers: [6],
        bondedRGroups: [2, 2, 2, 2, 2, 2],
      }
    },
    'two-separated-1': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T}|PEPTIDE2{[meI]}$$$$'},
      tgt: {
        simplePolymers: [5, 1],
        bondedRGroups: [1, 2, 2, 2, 1, 1],
      }
    },
    'two-separated-2': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz]}|PEPTIDE2{T.[meI]}$$$$'},
      tgt: {
        simplePolymers: [4, 2],
        bondedRGroups: [1, 2, 2, 1, 1, 1],
      }
    },
    'two-connected-1': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T}|PEPTIDE2{[meI]}$PEPTIDE1,PEPTIDE2,5:R2-1,R1$$$'},
      tgt: {
        simplePolymers: [5, 1],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    },
    'two-connected-2': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz]}|PEPTIDE2{T.[meI]}$PEPTIDE1,PEPTIDE2,4:R2-1,R1$$$'},
      tgt: {
        simplePolymers: [4, 2],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    }
  };

  for (const [testName, {src, tgt}] of Object.entries(tests)) {
    test(testName, async () => {
      const resHelm = new Helm(src.helm);

      const simplePolymers = resHelm.simplePolymers
        .map((sp) => sp.monomers.length);
      const totalMonomerCount = simplePolymers.reduce((a, b) => a + b, 0);
      expect(simplePolymers, tgt.simplePolymers);

      const bondedRGroups = wu.count(0).take(resHelm.bondedRGroupsMap.size)
        .map((i) => resHelm.bondedRGroupsMap.get(i)!.length).toArray();
      expect(totalMonomerCount, bondedRGroups.length);
      // expectArray(bondedRGroups, tgt.bondedRGroups);
    }, {skipReason: 'new tests'});
  }
});
