/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getHelmHelper, HelmInputBase, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerNumberingTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {HelmAtom, HelmMol} from '@datagrok-libraries/bio/src/helm/types';

import {_package} from '../package-test';

/**
 * Regression tests for the HELM-drawing interactivity used by the PolyTool
 * Enumeration dialog ({@link ../polytool/pt-enumerate-seq-dialog}). The dialog
 * lets the user pick positions by CLICKING monomers in the drawing and shows
 * substitute-monomer tooltips on HOVER, both routed through
 * `helmHelper.getHoveredAtom(e.offsetX, e.offsetY, mol, height)`.
 *
 * These broke after the hwe migration: the editor now paints an SVG with a
 * fit-to-surface transform, so (a) the view's `atom.p` was left in MODEL space
 * while the cursor reports SURFACE space, and (b) every SVG glyph is its own
 * hit target, so `offsetX/Y` was measured relative to the glyph under the
 * cursor instead of the host. Both are fixed in hwe (positions are baked to
 * surface space; the SVG is `pointer-events: none` so events fall through to
 * the host). The tests below click the ACTUAL rendered glyph (located via
 * `data-atom-id`, independent of `atom.p`) and assert the dialog's resolution
 * pipeline returns the monomer that was clicked.
 */

const W = 640;
const H = 220;
// A long-ish single chain so the fit transform scales the drawing UP to fill
// the pane (the maxScale: Infinity case) — i.e. surface coords differ from
// model coords and the monomers are spaced well apart.
const HELM = 'PEPTIDE1{A.C.D.E.F.G.H.I.K.L}$$$$V2.0';

category('PolyTool: Enumerate interaction', () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; // backup

  before(async () => {
    helmHelper = await getHelmHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  /** Builds the macromolecule input exactly as the enumeration dialog does,
   * mounts it at a fixed size, renders the HELM, and returns the handles a
   * caller needs to emulate pointer interaction. */
  function mountInput(): {input: HelmInputBase, host: HTMLElement, root: HTMLElement, mol: HelmMol} {
    const input = helmHelper.createHelmInput('Macromolecule', {
      editable: false,
      editorOptions: {drawOptions: {monomerNumbering: MonomerNumberingTypes.continuous}},
    });
    const root = ui.div([input.root], {
      style: {width: `${W}px`, height: `${H}px`, position: 'absolute', left: '0px', top: '0px'},
    });
    document.body.appendChild(root);
    const host = (input as any).viewer.host as HTMLElement;
    // Render the HELM, then pin a deterministic surface size (the input's own
    // ResizeObserver-driven sizing is async and would otherwise race the test).
    input.stringValue = HELM;
    (input as any).viewer.editor.setSize(W, H);
    return {input, host, root, mol: input.molValue};
  }

  /** Centre of a rendered atom glyph, in client (viewport) coordinates —
   * derived from the SVG element, NOT from `atom.p`, so it independently
   * validates the model→surface transform. */
  function glyphClientCentre(host: HTMLElement, atomId: number): {x: number, y: number} {
    const g = host.querySelector(`[data-atom-id="${atomId}"]`);
    expect(g != null, true, `glyph for atom ${atomId} must be rendered`);
    const shape = (g!.querySelector('.hw-atom__shape') ?? g!) as Element;
    const r = shape.getBoundingClientRect();
    return {x: r.left + r.width / 2, y: r.top + r.height / 2};
  }

  test('clicking a rendered monomer resolves to that monomer (all positions)', async () => {
    const {input, host, root, mol} = mountInput();
    try {
      expect(mol.atoms.length, 10);

      for (let k = 0; k < mol.atoms.length; k++) {
        const atomId = (mol.atoms[k].bio as any).id as number;
        const c = glyphClientCentre(host, atomId);

        // Dispatch on the host itself (the element `onClick` listens on). The
        // SVG is `pointer-events: none`, so a REAL click at this point lands on
        // the host too; dispatching here directly is robust to any shell
        // overlay that `elementFromPoint` might otherwise return. `offsetX/Y`
        // is then computed host-relative (= the baked surface position).
        let resolvedIdx = -2;
        let fired = false;
        const sub = input.onClick.subscribe((e: MouseEvent) => {
          fired = true;
          // EXACTLY the dialog's resolution (pt-enumerate-seq-dialog onClick).
          const a: HelmAtom | null = helmHelper.getHoveredAtom(e.offsetX, e.offsetY, mol, host.clientHeight);
          resolvedIdx = a ? mol.atoms.indexOf(a) : -1;
        });
        host.dispatchEvent(new MouseEvent('click', {clientX: c.x, clientY: c.y, bubbles: true}));
        sub.unsubscribe();

        expect(fired, true, `click on monomer ${k} must trigger the input's click event`);
        expect(resolvedIdx, k, `click on monomer ${k}'s glyph must resolve to position ${k}`);
      }
    } finally {
      (input as any).detach?.();
      root.remove();
    }
  });

  test('hovering a rendered monomer resolves to that monomer', async () => {
    const {input, host, root, mol} = mountInput();
    try {
      // Hover a middle monomer; the dialog shows the substitute tooltip for the
      // hovered position, so the hovered atom must resolve correctly.
      const k = 5;
      const atomId = (mol.atoms[k].bio as any).id as number;
      const c = glyphClientCentre(host, atomId);

      let resolvedIdx = -2;
      let fired = false;
      const sub = input.onMouseMove.subscribe((e: MouseEvent) => {
        fired = true;
        const a = helmHelper.getHoveredAtom(e.offsetX, e.offsetY, mol, host.clientHeight);
        resolvedIdx = a ? mol.atoms.indexOf(a) : -1;
      });
      host.dispatchEvent(new MouseEvent('mousemove', {clientX: c.x, clientY: c.y, bubbles: true}));
      sub.unsubscribe();

      expect(fired, true, 'mousemove over a monomer must trigger the input mousemove event');
      expect(resolvedIdx, k, 'hover over a monomer glyph must resolve to that position');
    } finally {
      (input as any).detach?.();
      root.remove();
    }
  });

  test('the editor SVG is transparent to pointer events', async () => {
    // The fix that makes `offsetX/Y` host-relative: with the SVG interactive,
    // a pointer over a monomer would target an SVG glyph and the dialog would
    // read glyph-relative coordinates (the reported bug). The deterministic
    // guarantee is the `pointer-events="none"` attribute on the SVG.
    const {input, host, root, mol} = mountInput();
    try {
      const svg = host.querySelector('svg');
      expect(svg != null, true, 'an SVG must be mounted');
      expect(svg!.getAttribute('pointer-events'), 'none');

      // Best-effort hit-test confirmation: only meaningful when our input is
      // actually the topmost element at that point (shell panels / the test
      // runner can sit above it in a live workspace). When the hit DOES land
      // inside our host, it must be the host — never an SVG glyph.
      const c = glyphClientCentre(host, (mol.atoms[3].bio as any).id as number);
      const top = document.elementFromPoint(c.x, c.y);
      if (top != null && host.contains(top))
        expect(svg!.contains(top) || top === svg, false, 'a hit inside the host must not be an SVG node');
    } finally {
      (input as any).detach?.();
      root.remove();
    }
  });

  test('atom positions are distinct and within the editor surface', async () => {
    const {input, root, mol} = mountInput();
    try {
      const seen = new Set<string>();
      for (const a of mol.atoms) {
        seen.add(`${Math.round(a.p.x)},${Math.round(a.p.y)}`);
        expect(a.p.x >= 0 && a.p.x <= W, true, `x in surface for ${a.p.x}`);
        expect(a.p.y >= 0 && a.p.y <= H, true, `y in surface for ${a.p.y}`);
      }
      // No collapse — every monomer occupies a distinct point (the degenerate
      // "all atoms at one spot" failure made hover/click unusable).
      expect(seen.size, mol.atoms.length);
    } finally {
      (input as any).detach?.();
      root.remove();
    }
  });
});
