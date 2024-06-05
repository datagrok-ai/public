import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HelmType} from '@datagrok/js-draw-lite/src/types/org';
import {OrgType} from '@datagrok/helm-web-editor/src/types/org-helm';
import {JSDraw2Module} from '../types';

declare const org: OrgType;
declare const JSDraw2: JSDraw2Module;

import {after, before, category, delay, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';
import {IMol} from '@datagrok/js-draw-lite/src/types/jsdraw2';

type TestTgtType = { atomCount: number, bondCount: number };

category('parseHelm', () => {
  const testData: { [testName: string]: { src: string, tgt: TestTgtType } } = {
    'PT': {
      src: 'PEPTIDE1{[L-hArg(Et,Et)].E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0',
      tgt: {atomCount: 6, bondCount: 5}
    },
    'RNA': {
      src: 'RNA1{[dhp]([m1A])p.r(T)p.r(T)p.r(G)p}$$$$V2.0',
      tgt: {atomCount: 12, bondCount: 11}
    }
  };

  for (const [testName, {src, tgt}] of Object.entries(testData)) {
    test(`Editor-${testName}`, async () => {
      _testParseHelmWithEditor(src, tgt);
    });
  }

  for (const [testName, {src, tgt}] of Object.entries(testData)) {
    test(`woDOM-${testName}`, async () => {
      _testParseHelmWithoutDOM(src, tgt);
    });
  }
});

function _testParseHelmWithEditor(src: string, tgt: TestTgtType): void {
  const io = org.helm.webeditor.IO;

  const editorHost = ui.div();
  const editor = new JSDraw2.Editor(editorHost, {viewonly: true});

  const plugin = new org.helm.webeditor.Plugin(editor);
  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, src, origin, null);

  const m: IMol<HelmType> = plugin.jsd.m;
  expect(m.atoms.length, tgt.atomCount);
  expect(m.bonds.length, tgt.bondCount);
}

function _testParseHelmWithoutDOM(src: string, tgt: TestTgtType): void {
  const io = org.helm.webeditor.IO;

  const molHandler = new JSDraw2.MolHandler();

  const plugin = new org.helm.webeditor.Plugin(molHandler);
  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, src, origin, null);

  const m: IMol<HelmType> = plugin.jsd.m;
  expect(m.atoms.length, tgt.atomCount);
  expect(m.bonds.length, tgt.bondCount);
}
