import {css, expect, flush, mount, need, row} from '../../lib/harness.mjs';

const COMPOUND = {id: 'V004044', name: 'Reference compound', smiles: 'CC(=O)Oc1ccccc1C(=O)O'};

async function build(mod) {
  const root = mount(need(mod, 'compoundSummary')({...COMPOUND}));
  await flush();
  const drawings = root.querySelectorAll('.eval-molecule');
  expect(drawings.length === 1, `the summary shows ${drawings.length} drawings of the molecule`);
  expect(drawings[0].dataset.smiles === COMPOUND.smiles, 'the drawing is not of the compound');
  return {root, drawing: drawings[0], structure: row('Structure', root)};
}

const px = (v) => /^([\d.]+)px$/.exec(v ?? '')?.[1];

export const checks = {
  'rows': async (mod) => {
    const {root} = await build(mod);
    for (const [caption, value] of [['ID', COMPOUND.id], ['Name', COMPOUND.name]]) {
      const r = row(caption, root);
      expect(r.text === value, `${caption} shows "${r.text}", not "${value}"`);
    }
  },
  'structure-in-value': async (mod) => {
    const {drawing, structure} = await build(mod);
    expect(structure.el.contains(drawing), 'the drawing is not in the Structure row');
    expect(drawing.closest('[data-u2-part="options"]') === null,
      'the drawing sits on the row\'s options rail (its side icons), not as its value');
  },
  'structure-has-room': async (mod) => {
    const {drawing, structure} = await build(mod);
    for (let el = drawing.parentElement; el && el !== structure.el.parentElement; el = el.parentElement) {
      for (const prop of ['height', 'max-height']) {
        const h = px(css.computed(el, prop));
        expect(h === undefined || +h >= 100, `the 100px drawing is held in a box ${h}px high (${el.className})`);
      }
    }
  },
};
