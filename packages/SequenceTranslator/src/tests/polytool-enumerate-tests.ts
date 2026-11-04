import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, awaitCheck, category, expect, test, expectArray} from '@datagrok-libraries/test/src/test';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes} from '../polytool/types';
import {doPolyToolEnumerateHelm} from '../polytool/pt-enumeration-helm';
import {polyToolEnumerateHelmUI} from '../polytool/pt-enumerate-seq-dialog';
import {PT_HELM_UI_DIALOG_ENUMERATION, PT_UI_HIGHLIGHT_MONOMERS} from '../polytool/const';

import {_package} from '../package-test';

category('PolyTool: Enumerate', () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; //backup

  before(async () => {
    helmHelper = await getHelmHelper(); // initialize JSDraw2 and org

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    // Clear settings to test default
    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const tests: {
    [testName: string]: { src: string, params: PolyToolEnumeratorParams, tgt: { seq: string, name: string }[] }
  } = {
    'single1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        placeholders: [
          {position: 4, monomers: ['K', 'P', 'F4COO']},
          {position: 6, monomers: ['Y', 'T']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5P'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5[F4COO]'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7T'},
      ]
    },
    'single-with-original': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        placeholders: [
          {position: 4, monomers: ['K', 'P', 'F4COO']},
          {position: 6, monomers: ['Y', 'T']},
        ],
        keepOriginal: true,
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: ''},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5P'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5[F4COO]'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7T'},
      ]
    },
    'parallel1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Parallel,
        placeholders: [
          {position: 1, monomers: ['D', 'L']},
          {position: 6, monomers: ['Y', 'T']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-[Tic]7T'},
      ],
    },
    'parallel-three-positions': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Parallel,
        placeholders: [
          {position: 1, monomers: ['D', 'L', 'K']},
          {position: 4, monomers: ['K', 'P', 'F4COO']},
          {position: 6, monomers: ['Y', 'T', 'A']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5K-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5P-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].K.W.G.[F4COO].L.A.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2K-P5[F4COO]-[Tic]7A'},
      ],
    },
    'parallel-with-original': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Parallel,
        placeholders: [
          {position: 1, monomers: ['D', 'L']},
          {position: 6, monomers: ['Y', 'T']},
        ],
        keepOriginal: true,
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: ''},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-[Tic]7T'},
      ],
    },
    'parallel-single-position': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Parallel,
        placeholders: [
          {position: 4, monomers: ['K', 'P', 'F4COO']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5P'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5[F4COO]'},
      ],
    },
    'matrix1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params:
        {
          type: PolyToolEnumeratorTypes.Matrix,
          placeholders: [
            {position: 1, monomers: ['D', 'L']},
            {position: 4, monomers: ['K', 'P', 'F4COO']},
            {position: 6, monomers: ['Y', 'T']},
          ]
        },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5K-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5K-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5P-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5P-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5[F4COO]-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5[F4COO]-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5K-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5K-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5P-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5P-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5[F4COO]-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5[F4COO]-[Tic]7T'},
      ],
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const res = doPolyToolEnumerateHelm(testData.src, '', testData.params);
      expectArray(res, testData.tgt.map((r) => [r.seq, r.name]));
    });
  }

  // The box was hidden for macromolecules, where its value (restored from the dialog's history) still applied, so
  // every enumerated molecule could come out highlighted with no box to untick
  test('dialog: Highlight monomers shows with To atomic level for a macromolecule', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('seq', ['PEPTIDE1{A.G.C.D}$$$$'])]);
    await grok.data.detectSemanticTypes(df);
    await polyToolEnumerateHelmUI(df.cell(0, 'seq'));
    const dlg = DG.Dialog.getOpenDialogs().find((d) => d.title === PT_HELM_UI_DIALOG_ENUMERATION);
    expect(!!dlg, true, 'the dialog did not open');
    try {
      const inputRoot = (caption: string): HTMLElement => Array.from(dlg!.root.querySelectorAll('div.ui-input-root'))
        .find((r) => r.querySelector('label.ui-input-label')?.textContent?.trim() === caption) as HTMLElement;
      const highlight = inputRoot(PT_UI_HIGHLIGHT_MONOMERS);
      const toAtomicLevel = inputRoot('To atomic level').querySelector('input[type="checkbox"]') as HTMLInputElement;
      // the input reports a change after the click's event, not within it
      if (toAtomicLevel.checked)
        toAtomicLevel.click();
      await awaitCheck(() => highlight.style.display === 'none', 'shown with To atomic level off', 3000);
      toAtomicLevel.click();
      expect(toAtomicLevel.checked, true, 'To atomic level did not switch on');
      await awaitCheck(() => highlight.style.display === '', 'hidden with To atomic level on', 3000);
    } finally {
      dlg!.close();
    }
  });

  // bio-edge 2026-11-02, pinned (not fixed): a placeholder past the sequence's end, or negative, is ignored and each
  // of its monomers gives a copy of the source (Enumerate Single HELM Sequence takes positions as they come).
  // Listed: an error or no rows would change what its callers get.
  test('position-out-of-range-pinned', async () => {
    const src = 'PEPTIDE1{A.C.D}$$$$V2.0';
    for (const position of [10, -1]) {
      const res = doPolyToolEnumerateHelm(src, '', {
        type: PolyToolEnumeratorTypes.Single, placeholders: [{position, monomers: ['G', 'K']}],
      });
      expectArray(res.map(([seq]) => seq), [src, src]);
    }
  });
});
