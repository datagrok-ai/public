import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {splitterAsHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {ISeqGraphInfoHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

category('helm', () => {
  const tests = {
    'single-linear': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T.[meI]}$$$$'},
      tgt: {
        simplePolymers: [6], connections: [],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    },
    'single-cyclized-C-2-2': {
      src: {helm: 'PEPTIDE1{R.F.C.Y.G.H.[GGaz].C.T.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-8:R3$$$'},
      tgt: {
        simplePolymers: [10], connections: [[['PEPTIDE1', 3, 'R3'], ['PEPTIDE1', 8, 'R3']]],
        bondedRGroups: [1, 2, 3, 2, 2, 2, 2, 3, 2, 1],
      }
    },
    'single-cyclized-C-1-1': {
      src: {helm: 'PEPTIDE1{F.C.Y.G.H.[GGaz].C.[meI]}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$'},
      tgt: {
        simplePolymers: [8], connections: [[['PEPTIDE1', 2, 'R3'], ['PEPTIDE1', 7, 'R3']]],
        bondedRGroups: [1, 3, 2, 2, 2, 1, 3, 1],
      }
    },
    'single-cyclized-C-0-0': {
      src: {helm: 'PEPTIDE1{C.Y.G.H.[GGaz].C}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$'},
      tgt: {
        simplePolymers: [6], connections: [[['PEPTIDE1', 1, 'R3'], ['PEPTIDE1', 6, 'R3']]],
        bondedRGroups: [2, 2, 2, 2, 2, 2],
      }
    },
    'two-separated-5-1': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T}|PEPTIDE2{[meI]}$$$$'},
      tgt: {
        simplePolymers: [5, 1], connections: [],
        bondedRGroups: [1, 2, 2, 2, 1, 0],
      }
    },
    'two-separated-1-5': {
      src: {helm: 'PEPTIDE1{[meI]}|PEPTIDE2{R.F.Y.[GGaz].T}$$$$'},
      tgt: {
        simplePolymers: [1, 5], connections: [],
        bondedRGroups: [0, 1, 2, 2, 2, 1],
      }
    },
    'two-separated-4-2': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz]}|PEPTIDE2{T.[meI]}$$$$'},
      tgt: {
        simplePolymers: [4, 2], connections: [],
        bondedRGroups: [1, 2, 2, 1, 1, 1],
      }
    },
    'two-connected-1': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz].T}|PEPTIDE2{[meI]}$PEPTIDE1,PEPTIDE2,5:R2-1:R1$$$'},
      tgt: {
        simplePolymers: [5, 1], connections: [[['PEPTIDE1', 5, 'R2'], ['PEPTIDE2', 1, 'R1']]],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    },
    'two-connected-2': {
      src: {helm: 'PEPTIDE1{R.F.Y.[GGaz]}|PEPTIDE2{T.[meI]}$PEPTIDE1,PEPTIDE2,4:R2-1:R1$$$'},
      tgt: {
        simplePolymers: [4, 2], connections: [[['PEPTIDE1', 4, 'R2'], ['PEPTIDE2', 1, 'R1']]],
        bondedRGroups: [1, 2, 2, 2, 2, 1],
      }
    },
    'two-cyclized-1-9': {
      src: {helm: 'PEPTIDE1{[meI]}|PEPTIDE2{R.F.[GGaz].T.G.H.F.Y.P}$PEPTIDE2,PEPTIDE2,3:R3-9:R2|PEPTIDE2,PEPTIDE1,3:R4-1:R1$$$V2.0'},
      tgt: {
        simplePolymers: [1, 9],
        connections: [
          [['PEPTIDE2', 3, 'R3'], ['PEPTIDE2', 9, 'R2']],
          [['PEPTIDE2', 3, 'R4'], ['PEPTIDE1', 1, 'R1']]],
        bondedRGroups: [1, 1, 2, 4, 2, 2, 2, 2, 2, 1],
      }
    }
  };

  // The converter's own HELM parser went with the JavaScript assembly (Bio on Crux): the same cases now check the
  // splitter (hwe adapter) whose chains and connections feed the assembler
  for (const [testName, {src, tgt}] of Object.entries(tests)) {
    test(testName, async () => {
      const seqSS = splitterAsHelm(src.helm);
      const gi = seqSS.graphInfo as ISeqGraphInfoHelm;
      const starts = gi.disjointSeqStarts;

      const simplePolymers = wu.count(0).take(starts.length)
        .map((i) => (i + 1 < starts.length ? starts[i + 1] : seqSS.length) - starts[i]).toArray();
      const totalMonomerCount = simplePolymers.reduce((a, b) => a + b, 0);
      expectArray(simplePolymers, tgt.simplePolymers);
      expect(totalMonomerCount, seqSS.length);

      const chainIds = gi.chainIds ?? [];
      const connections = gi.connections
        .map((c) => {
          return [
            [chainIds[c.seqIndex1], c.monomerIndex1 + 1, `R${c.rGroup1}`],
            [chainIds[c.seqIndex2], c.monomerIndex2 + 1, `R${c.rGroup2}`]];
        });

      expectArray(connections, tgt.connections);
    });
  }
});
