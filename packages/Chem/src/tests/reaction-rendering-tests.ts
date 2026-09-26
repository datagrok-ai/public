import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {before, category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {prepareReactionDepiction} from '../rendering/reaction-depiction';
import {BRANCH_DELIMITER, RDKitReactionRenderer} from '../rendering/rdkit-reaction-renderer';

category('Reaction rendering', () => {
  let rdkit: RDModule;
  before(async () => {rdkit = await grok.functions.call('Chem:getRdKitModule');});

  function expectAligned(molblocks: Map<string, string>, scaffold: string, product: string): void {
    const reference = rdkit.get_mol(scaffold);
    const target = rdkit.get_mol(product);
    const a = rdkit.get_mol(molblocks.get(reference.get_smiles())!);
    const b = rdkit.get_mol(molblocks.get(target.get_smiles())!);
    try {
      const match = JSON.parse(b.get_substruct_match(a)).atoms;
      expect(!!match, true);
      const ac = JSON.parse(a.get_json()).molecules[0].conformers[0].coords;
      const bc = JSON.parse(b.get_json()).molecules[0].conformers[0].coords;
      for (let i = 1; i < match.length; i++) {
        for (let axis = 0; axis < 2; axis++)
          expectFloat(ac[i][axis] - ac[0][axis], bc[match[i]][axis] - bc[match[0]][axis], 0.001);
      }
    } finally {
      reference.delete(); target.delete(); a.delete(); b.delete();
    }
  }

  test('esterification aligns the retained acid instead of water or ethanol', async () => {
    const molblocks = new Map<string, string>();
    const depiction = prepareReactionDepiction(rdkit, 'CC(=O)O.CCO>>CCOC(C)=O.O', molblocks);
    expect(depiction !== null, true);
    expectAligned(molblocks, 'CC(=O)O', 'CCOC(C)=O');
    const reaction = rdkit.get_rxn(depiction!);
    try {expect(!!reaction, true);} finally {reaction?.delete();}
  });

  test('aspirin keeps its scaffold and repeated intermediate across steps', async () => {
    const acid = 'O=C(O)c1ccccc1O';
    const aspirin = 'CC(=O)Oc1ccccc1C(=O)O';
    const molblocks = new Map<string, string>();
    expect(prepareReactionDepiction(rdkit, `${acid}.CC(=O)OC(C)=O>>${aspirin}`, molblocks) !== null, true);
    expectAligned(molblocks, acid, aspirin);
    const first = new Map(molblocks);
    expect(prepareReactionDepiction(rdkit, `${aspirin}>>${acid}.CC(=O)O`, molblocks) !== null, true);
    expectAligned(molblocks, acid, aspirin);
    for (const [key, block] of first)
      expect(molblocks.get(key), block);
  });

  test('preserves stereochemistry, charges and atom maps in cached molblocks', async () => {
    const substrate = '[NH3+:1][C@@H:2](C)C(=O)[O-:3]';
    const molblocks = new Map<string, string>();
    expect(prepareReactionDepiction(rdkit, `${substrate}>>${substrate}`, molblocks) !== null, true);
    for (const [smiles, block] of molblocks) {
      const molecule = rdkit.get_mol(block);
      try {expect(molecule.get_smiles(), smiles);} finally {molecule.delete();}
    }
  });

  test('keeps query SMARTS and supplied CX coordinates on their original rendering path', async () => {
    const molblocks = new Map<string, string>();
    for (const reaction of ['[C,N:1]-[O;H1:2]>>[C,N:1]=[O:2]', 'CC>>CO |(0,0,;1,0,;0,0,;1,0,)|'])
      expect(prepareReactionDepiction(rdkit, reaction, molblocks), null);
    expect(molblocks.size, 0);
  });

  test('accepts molecule blocks without losing their header or supplied coordinates', async () => {
    const mol = rdkit.get_mol('CC(=O)O');
    try {
      mol.set_new_coords(false);
      const block = mol.get_molblock();
      const molblocks = new Map<string, string>();
      const input = `${block}>>${block}`;
      expect(prepareReactionDepiction(rdkit, input, molblocks) !== null, true);
      expect(molblocks.get(mol.get_smiles()), block);
      const renderer = new RDKitReactionRenderer(rdkit);
      expect(renderer.renderToCanvas(ui.canvas(600, 150), input), true);
    } finally {
      mol.delete();
    }
  });

  test('keeps agents separate from reactants and products', async () => {
    const molblocks = new Map<string, string>();
    const depiction = prepareReactionDepiction(rdkit, 'CC(=O)O.CCO>O>CCOC(C)=O', molblocks);
    expect(depiction !== null, true);
    const rxn = rdkit.get_rxn(depiction!);
    try {
      expect(rxn.get_svg(600, 150).includes('<svg'), true);
      expectAligned(molblocks, 'CC(=O)O', 'CCOC(C)=O');
    } finally {
      rxn.delete();
    }
  });

  test('reuses reaction images on repaint and coordinates on resize', async () => {
    const renderer = new RDKitReactionRenderer(rdkit);
    const canvas = ui.canvas(900, 240);
    const route = 'CC(=O)O.CCO>>CCOC(C)=O>>CC(=O)O';
    expect(renderer.renderToCanvas(canvas, route), true);
    const renders = renderer.canvasCounter;
    expect(renderer.renderToCanvas(canvas, route), true);
    expect(renderer.canvasCounter, renders);
    expect(renderer.renderToCanvas(canvas, route, 450, 240), true);
    expect(renderer.canvasCounter > renders, true);
  });

  test('wrapped steps stay inside their panel at high DPI and preserve the canvas transform', async () => {
    const renderer = new RDKitReactionRenderer(rdkit);
    const canvas = ui.canvas(1500, 800);
    const ctx = canvas.getContext('2d')!;
    ctx.scale(2, 2);
    const steps = ['CCO>>CC=O', 'CC=O>>CC(=O)O', 'N>>CN', 'CN.CC(=O)O>>CNC(C)=O'];
    renderer._drawMultiStepReaction(80, 60, 1300, 600, canvas, steps.map((s) => [s]), 2);
    expect(ctx.getTransform().a, 2);
    const pixels = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
    let ink = 0;
    for (let y = 0; y < canvas.height; y++) {
      for (let x = 0; x < canvas.width; x++) {
        const alpha = pixels[(y * canvas.width + x) * 4 + 3];
        if (x < 80 || x >= 1380 || y < 60 || y >= 660)
          expect(alpha, 0);
        else
          ink += alpha > 0 ? 1 : 0;
      }
    }
    expect(ink > 100, true);
    expect(renderer.renderToCanvas(canvas, steps.join(BRANCH_DELIMITER), 5, 5), true);
  });
});
