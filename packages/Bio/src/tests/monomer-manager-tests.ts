import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, category, expect, expectArray, expectObject, test} from '@datagrok-libraries/test/src/test';
import {Monomer, RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {Control, Scope} from '@datagrok-libraries/u2';

import {CapGroups, COMMON_CAP_GROUPS} from '../utils/monomer-lib/monomer-manager/cap-groups';
import {RGroupsEditor} from '../utils/monomer-lib/monomer-manager/rgroups-editor';
import {MonomerEditor} from '../utils/monomer-lib/monomer-manager/monomer-editor';
import {MonomerManager} from '../utils/monomer-lib/monomer-manager/monomer-manager';
import {MonomerLibManager} from '../utils/monomer-lib/lib-manager';

const rg = (label: string, smiles: string, name: string, alternateId: string): RGroup =>
  ({label, capGroupSMILES: smiles, capGroupName: name, alternateId});

const canonicalize = (smiles: string): string | null =>
  grok.chem.checkSmiles(smiles) ? grok.chem.convert(smiles, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) : null;

/** Every way the production libraries write a cap (HELM core, Regmol, monomersforDataGrok), and the
 * broken forms the manager's fixers have met. */
const MESSY: [RGroup, string | null, boolean][] = [
  [rg('R1', '[*:1][H]', 'H', 'R1-H'), 'H', true],
  [rg('R1', '[H][*:1]', 'H', 'R1-H'), 'H', true],
  [rg('R2', 'O[*:2]', 'OH', 'R2-OH'), 'OH', true],
  [rg('R2', '[O:2]', '', ''), 'OH', true],
  [rg('R2', 'O[R2]', '', ''), 'OH', true],
  [rg('R2', '[OH][*:2]', '', ''), 'OH', true],
  [rg('R3', '[*:3]', 'Br', 'R3-Br'), 'Br', true],
  [rg('R3', '[*:3]Azide', 'Azide', 'R3-Azide'), 'Azide', true],
  [rg('R3', '[*:3]Allyl', 'Allyl', 'R3-Allyl'), 'Allyl', true],
  [rg('R3', 'C=C[*:3]', 'Allyl', 'R3-Allyl'), 'Allyl', true],
  [rg('R1', 'C[*:1]=C', 'Allyl', 'R1-Allyl'), 'Allyl', true],
  [rg('R1', 'C[*:1]#C', 'Alkyne', 'R1-Alkyne'), 'Alkyne', true],
  [rg('R3', '[*:3]Alkynyl', 'Alkynyl', 'R3-Alkynyl'), 'Alkyne', true],
  [rg('R2', 'Cl[*:2]', 'Cl', 'R2-Cl'), 'Cl', true],
  [rg('R3', 'I[*:3]', 'I', 'R3-I'), 'I', true],
  [rg('R1', '', 'hydrogen', ''), 'H', true],
  [rg('R1', '', '', 'R1-H'), 'H', true],
  [rg('R1', 'C(C)(=O)[*:1]', '', ''), 'Ac', true],
  [rg('R1', 'X[*:1]', 'X', 'R1-X'), 'X', false],
  [rg('R7', 'CC(C)(C)OC(=O)[*:7]', 'Boc', 'R7-Boc'), 'Boc', false],
  [rg('R5', '', '', ''), null, false],
];

category('monomerManager: cap groups', () => {
  test('infers the cap of every way a library writes it', async () => {
    for (const [def, name, known] of MESSY) {
      const inferred = CapGroups.infer(def, canonicalize);
      const what = JSON.stringify(def);
      expect(inferred?.cap.name ?? null, name, `${what}: cap`);
      if (inferred !== null)
        expect(inferred.known, known, `${what}: known`);
    }
  });

  test('a name over different chemistry keeps the chemistry', async () => {
    const inferred = CapGroups.infer(rg('R3', 'C=CC[*:3]', 'Allyl', 'R3-Allyl'), canonicalize)!;
    expect(inferred.known, false);
    expect(inferred.cap.smiles, 'C=CC[*]');
    expect(inferred.note, 'not the usual Allyl');
    const conflict = CapGroups.infer(rg('R8', 'Cl[*:8]', 'Br', 'R8-Br'), canonicalize)!;
    expect(conflict.cap.name, 'Cl');
    expect(conflict.note, 'the name says Br');
    // a typed SMILES that is no chemistry stays what will be saved, flagged, instead of the name beside it
    const typo = CapGroups.infer(rg('R3', 'C1CC[*:3', 'H', 'R3-H'), canonicalize)!;
    expect(typo.known, false);
    expect(typo.note, 'the SMILES is not valid');
  });

  test('capping joins each cap at the atom it closes, whatever its size', async () => {
    const ala = '[*:1]N[C@@H](C)C(=O)[*:2]';
    const same = (a: string, b: string) => expect(canonicalize(a), canonicalize(b), `${a} vs ${b}`);
    same(CapGroups.capped(ala, [rg('R1', '[*:1][H]', 'H', 'R1-H'), rg('R2', 'O[*:2]', 'OH', 'R2-OH')],
      canonicalize), 'N[C@@H](C)C(=O)O');
    same(CapGroups.capped(ala, [rg('R1', 'CC(=O)[*:1]', 'Ac', 'R1-Ac'), rg('R2', 'CO[*:2]', 'OMe', 'R2-OMe')],
      canonicalize), 'CC(=O)N[C@@H](C)C(=O)OC');
    // nothing defined, or a placeholder: closed with H
    same(CapGroups.capped(ala, [rg('R2', 'X[*:2]', 'X', 'R2-X')], canonicalize), 'N[C@@H](C)C=O');
  });

  test('a known cap makes the R group the libraries write', async () => {
    expectObject(CapGroups.rgroup(CapGroups.byName('OH')!, 2),
      {capGroupSMILES: 'O[*:2]', alternateId: 'R2-OH', capGroupName: 'OH', label: 'R2'});
    expectObject(CapGroups.rgroup(CapGroups.byName('H')!, 1),
      {capGroupSMILES: '[*:1][H]', alternateId: 'R1-H', capGroupName: 'H', label: 'R1'});
    expect(COMMON_CAP_GROUPS.length, 15);
    expectArray(COMMON_CAP_GROUPS.slice(0, 3).map((c) => c.name), ['H', 'OH', 'Allyl']);
  });
});

category('monomerManager: R-groups editor', () => {
  let live = 0;
  const editors: RGroupsEditor[] = [];
  const make = (): RGroupsEditor => {
    live = Scope.liveCount;
    const editor = new RGroupsEditor({canonicalize, isValidSmiles: (s) => grok.chem.checkSmiles(s)});
    document.body.append(editor.root);
    editors.push(editor);
    return editor;
  };
  const dispose = (editor: RGroupsEditor) => {
    editor.root.remove();
    editor.dispose();
    editors.splice(editors.indexOf(editor), 1);
    expect(Scope.liveCount, live, 'scopes left behind by the editor');
  };
  const row = (editor: RGroupsEditor, label: string) =>
    editor.root.querySelector<HTMLElement>(`.monomer-rgroup-row[data-u2-item="${label}"]`)!;

  after(async () => {
    for (const editor of editors.splice(0)) {
      editor.root.remove();
      editor.dispose();
    }
  });

  test('saves common caps canonically and keeps the others as defined', async () => {
    const editor = make();
    editor.setDefinitions([rg('R1', '[H][*:1]', 'H', 'R1-H'), rg('R2', '[O:2]', '', ''),
      rg('R3', 'X[*:3]', 'X', 'R3-X'), rg('R4', '', '', '')]);
    const saved = editor.rgroups.value;
    expectObject(saved[0], {capGroupSMILES: '[*:1][H]', capGroupName: 'H', alternateId: 'R1-H', label: 'R1'});
    expectObject(saved[1], {capGroupSMILES: 'O[*:2]', capGroupName: 'OH', alternateId: 'R2-OH', label: 'R2'});
    expectObject(saved[2], {capGroupSMILES: 'X[*:3]', capGroupName: 'X', label: 'R3'});
    expect(editor.problem.value, 'R3: X is a placeholder, not a cap group — pick the cap R3 gets (usually H)');
    const select = row(editor, 'R3').querySelector('select')!;
    expect(select.value, '~own', 'the cap the library names is the chosen one');
    expect(select.options[select.options.length - 1].textContent, 'X', 'and is appended to the list');
    expect(row(editor, 'R4').querySelector('select')!.value, '', 'nothing to infer: the empty choice');
    dispose(editor);
  });

  test('a pick and the four fields both edit the R group', async () => {
    const editor = make();
    editor.setDefinitions([rg('R1', '[*:1][H]', 'H', 'R1-H'), rg('R2', '', '', '')]);
    const select = row(editor, 'R2').querySelector('select')!;
    select.value = 'Br';
    select.dispatchEvent(new Event('change'));
    expectObject(editor.rgroups.value[1], {capGroupSMILES: 'Br[*:2]', capGroupName: 'Br', alternateId: 'R2-Br'});
    expect(editor.problem.value, null);
    expect(row(editor, 'R2').querySelector('.monomer-rgroup-info')!.textContent, 'bromo · Br[*:2]');

    const fields = row(editor, 'R2').querySelector<HTMLElement>('.monomer-rgroup-fields')!;
    expect(fields.hidden, true, 'the fields wait for the pencil');
    row(editor, 'R2').querySelector<HTMLButtonElement>('[data-u2="icon-button"]')!.click();
    expect(fields.hidden, false);
    const [smiles, name, alternateId] = Array.from(fields.querySelectorAll('input'));
    for (const [input, value] of [[smiles, '[*:2]C(F)(F)F'], [name, 'CF3'], [alternateId, 'R2-CF3']] as const) {
      input.value = value;
      input.dispatchEvent(new Event('change'));
    }
    expectObject(editor.rgroups.value[1],
      {capGroupSMILES: '[*:2]C(F)(F)F', capGroupName: 'CF3', alternateId: 'R2-CF3'});
    expect(select.value, '~own');
    expect(select.options[select.options.length - 1].textContent, 'CF3');
    dispose(editor);
  });

  test('follows the structure\'s labels and reads caps set by name', async () => {
    const editor = make();
    editor.setDefinitions([rg('R1', '[*:1][H]', 'H', 'R1-H'), rg('R2', 'Cl[*:2]', 'Cl', 'R2-Cl')]);
    editor.setLabels(['R1', 'R2', 'R3']);
    expectArray(editor.rgroups.value.map((r) => r.capGroupName), ['H', 'Cl', 'H']);
    editor.setLabels(['R2']);
    expectArray(editor.labels, ['R2']);
    expect(editor.root.querySelectorAll('.monomer-rgroup-row').length, 1);
    editor.setCap('R2', 'OMe');
    expectObject(editor.rgroups.value[0], {capGroupSMILES: 'CO[*:2]', capGroupName: 'OMe'});
    let thrown = '';
    try {
      editor.setCap('R2', 'nonsense');
    } catch (e) {
      thrown = (e as Error).message;
    }
    expect(thrown.startsWith('"nonsense" is not one of the common caps'), true, thrown);
    dispose(editor);
  });
});

category('monomerManager: editor', () => {
  const ALA: Monomer = {
    symbol: 'bddAla', name: 'test alanine', molfile: '', smiles: '[*:1]N[C@@H](C)C([*:2])=O', polymerType: 'PEPTIDE',
    monomerType: 'Backbone', id: 7, naturalAnalog: 'A', author: '', createDate: null,
    rgroups: [rg('R1', '[*:1][H]', 'H', 'R1-H'), rg('R2', 'O[*:2]', 'OH', 'R2-OH')],
    meta: {note: 'x', mw: 89.09, d: true},
  };

  test('holds a monomer, says why it cannot be saved and goes without a trace', async () => {
    const libManager = await MonomerLibManager.getInstance();
    const live = Scope.liveCount;
    const editor = new MonomerEditor({libManager, activeLib: () => null, monomersDf: () => undefined,
      refreshTable: async () => {}, openLibrary: async () => {}, createLibrary: () => {}}, ['a.json'], 'a.json');
    document.body.append(editor.root);
    try {
      editor.setMonomer(ALA);
      expectObject(editor.getFields(), {polymerType: 'PEPTIDE', monomerType: 'Backbone', symbol: 'bddAla',
        name: 'test alanine', id: 7, naturalAnalog: 'A'});
      expect(editor.problem.value, null);
      expectObject(editor.meta.value.value, {note: 'x', mw: '89.09', d: 'true'});
      const saved = editor.buildMonomer().meta!;
      expect(saved.mw, 89.09, 'a meta value left as it was keeps its type');
      expect(saved.d, true);
      expect('colors' in saved, false, 'colors nobody changed are not written');
      editor.setFields({name: ''});
      expect(editor.problem.value, 'Monomer Name is required');
      editor.setFields({name: 'test alanine'});
      editor.rgroups.setCap('R2', 'NH2');
      expectArray(editor.buildMonomer().rgroups.map((r) => r.capGroupName), ['H', 'NH2']);
      editor.setEmptyMonomer();
      expect(editor.problem.value, 'Draw the monomer\'s structure');
      expect(editor.id.value.value, 0, 'a new monomer starts at ID 0, as the view does');
      const save = editor.root.querySelector<HTMLButtonElement>('.monomer-manager-save button')!;
      expect(save.disabled, true, 'Save waits for a structure');
    } finally {
      editor.root.remove();
      editor.dispose();
    }
    expect(Scope.liveCount, live, 'scopes left behind by the editor');
  });

  test('the Manage Monomers view edits through the copilot fields and disposes its editor on close', async () => {
    // the manager of this bundle — the package's own bundle keeps a singleton of its own
    const manager = await MonomerManager.getInstance();
    const view = grok.shell.v = await manager.getViewRoot();
    await awaitCheck(() => view.root.querySelector('.monomer-manager-editor') !== null,
      'the Manage Monomers view did not dock its editor', 30000);
    try {
      const editor = Control.forElement(view.root.querySelector('.monomer-manager-editor')) as MonomerEditor;
      let disposed = false;
      editor.own(() => disposed = true);
      const unrecognized = manager.setMonomerFormFields({molecule: ALA.smiles, symbol: 'bddAla', name: 'test alanine',
        polymerType: 'PEPTIDE', R2: 'NH2', color: 'red'});
      expectArray(unrecognized, ['color']);
      expectArray(editor.rgroups.rgroups.value.map((r) => r.capGroupName), ['H', 'NH2']);
      expect(manager.validateMonomerForm(), null);
      view.close();
      await awaitCheck(() => disposed, 'the editor outlived its view', 10000);
    } finally {
      Array.from(grok.shell.tableViews).find((v) => v.name === MonomerManager.VIEW_NAME)?.close();
    }
  });
});
