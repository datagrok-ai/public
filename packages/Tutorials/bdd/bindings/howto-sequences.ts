/* The Sequences how-to clips (features/howto/09-sequences): Bio's sample tables and the few Bio steps they use,
   copied from Bio's own bdd (bindings/elements.ts, steps.ts, demo-atomic.ts, demo-search.ts), which the Tutorials
   project does not load. Everything else is library vocabulary. */
import type {Page} from '@playwright/test';
import {Given, Then, When, dataset, element, kind} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {expect, expectCustomEvent, pollMs, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

dataset('FASTA_PT', {path: 'System:AppData/Bio/samples/FASTA_PT.csv', aliases: ['fasta-pt'],
  description: '100 peptides of one family, 31 to 38 residues, unaligned: id, sequence (fasta); MTMFC near the start and EAM.K near the end at shifting offsets'});
dataset('HELM_sample', {path: 'System:AppData/Bio/samples/HELM.csv', aliases: ['helm-sample'],
  description: '540 peptides: HELM (helm notation, multichar monomers), Activity'});
dataset('FASTA_sample', {path: 'System:AppData/Bio/samples/FASTA.csv', aliases: ['fasta-sample'],
  description: '64 UniProt peptides: Entry, Length, UniProtKB, Sequence (fasta), Activity, Cluster'});
dataset('MSA_sample', {path: 'System:AppData/Bio/samples/MSA.csv', aliases: ['msa-sample'],
  description: '540 aligned peptides: MSA (separator "/", multichar monomers), Activity'});

/* Antibody numbering (Bio > Annotate > Apply Numbering Scheme...), as Bio's demo bindings have it (demo-antibody.ts). */
dataset('antibodies', {path: 'System:AppData/Bio/samples/antibodies.csv',
  description: '493 antibodies: AntibodyHC and AntibodyLC (fasta peptides ~130–215 aa), Antigen, Y'});

/** Bio's header sets a new long sequence column 400 px wide shortly after it first paints the column's header, over a
 * width set meanwhile: the width is set once the header shows and set again whenever it changes in the next 400 ms. */
export const columnStaysWide = When('the {string} column of the grid stays {int} px wide',
  async (page: Page, column: string, width: number) => {
    silent(page);
    await page.evaluate(async ([c, w]) => {
      const col = grok.shell.tv?.grid?.columns.byName(c);
      if (!col)
        throw new Error(`the grid of the current view has no "${c}" column`);
      const until = Date.now() + 400;
      while (Date.now() < until) {
        if (col.width !== w)
          col.width = w;
        await new Promise((r) => requestAnimationFrame(r));
      }
    }, [column, width] as [string, number]);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a clip leaves out: the column\'s width, held past the header\'s own 400 px default'});

/** The region band (Bio's annotation track) is the lowest track of the header, right under the WebLogo whose
 * `position <name> of <column> header` areas the grid reports: the pointer goes into the band under that position. */
export const hoverRegionBand = When('user hovers over the region track under the {string} area of {widget}',
  async (page: Page, area: string, target: ElementRef) => {
    const logo = await viewers.hitArea(page, target, area, true);
    const x = logo.x + logo.width / 2;
    const y = logo.y + logo.height + 4 + 7;
    await page.mouse.move(x - 3, y - 3);
    await page.mouse.move(x, y);
    await expect(page.locator('.d4-tooltip').filter({visible: true}), `the tooltip of the region under ${area}`).toHaveCount(1);
  }, {tier: 'ui', description: 'the annotation track 4 px under the WebLogo stack of that position; done when its tooltip shows'});

/** Diversity Search picks its subset from a random start, so the rows it lists vary from run to run; how many
 * groups of a column they span is what holds. */
export const subsetSpans = Then('the peptides {widget} lists should come from {int} different {string} values',
  async (page: Page, target: ElementRef, count: number, column: string) => {
    await expect.poll(async () => {
      const subset = String(await viewers.readingOf(page, target, 'subset'));
      return page.evaluate(([ids, c]) => {
        const col = grok.shell.t?.col(c);
        if (!col)
          return `the table has no "${c}" column`;
        return new Set(ids.split(',').filter((s) => s !== '').map((i) => String(col.get(Number(i))))).size;
      }, [subset, column] as [string, string]);
    }, {timeout: pollMs(10000), message: `the "${column}" values of the rows ${target.phrase} lists`}).toBe(count);
  }, {description: 'the viewer\'s "subset" reading (row indexes) looked up in that column of the current table'});

/* The Peptides package's SAR (Bio > Analyze > SAR...), as Peptides' bdd has it (bindings/elements.ts, steps.ts). */
dataset('peptides', {path: 'System:DemoFiles/bio/peptides.csv',
  description: '647 aligned separator peptides, ID and IC50; 17 positions including NH2 and COOH'});

export const peptidesReady = Given('the Peptides package is initialized', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => { await grok.functions.call('Peptides:initPeptides'); });
}, {tier: 'api', description: 'waits for the sequence helper and monomer library initialization'});

export const sarAnalysisReady = Then('the SAR analysis should be ready', async (page: Page) => {
  const event = await expectCustomEvent(page, 'peptides-sar-ready', 180000) as {table: string};
  expect(event.table).toBe(await page.evaluate(() => grok.shell.t.name));
}, {description: 'consumes peptides-sar-ready after launch, including MCL completion'});

/* The monomer manager's editor (Bio > Manage > Monomers), as Bio's bdd names it (bindings/monomer-form.ts). */
element('monomer form', {selector: '[data-u2-name="monomerEditor"]',
  description: 'the editor of the Manage Monomers view: the library, the sketcher, the monomer\'s symbol, name and types, a chip per R group'});
kind('R-group chip', {selector: '.monomer-rgroup-chip', match: ['item'],
  description: 'an R group of the monomer in the editor, by its label ("R2"), reading "R2 OH"; a click opens its menu of the common caps'});

/** The monomer editor's sketcher has its toolbar (Ketcher), its canvas, or is initialized (any other backend). */
export const monomerEditorReady = Then('the monomer sketcher of the Manage Monomers view should be ready', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const root = document.querySelector('.monomer-manager-sketcher');
    if (!root)
      return 'no sketcher';
    const ketcher = root.querySelector('.Ketcher-root');
    if (ketcher)
      return ketcher.querySelectorAll('button').length > 5 ? 'ready' : 'mounting';
    const shown = root.querySelector('.chem-sketcher-host > *');
    const sketcher = shown ? (window as any).DG.Widget.find(shown) : null;
    return root.querySelector('canvas') || sketcher?.isInitialized ? 'ready' : 'mounting';
  }), {timeout: pollMs(150000), message: 'the monomer editor sketcher'}).toBe('ready');
}, {description: 'the sketcher inside .monomer-manager-sketcher has its toolbar, its canvas, or is initialized'});

/** Bio initializes on its first call (RDKit, the monomer libraries, the sequence helper), and the platform holds
 * every call of the package until then: a service getter awaits it, so a clip's first command does not pay for it. */
export const bioReady = Given('the Bio package is initialized', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => { await grok.functions.call('Bio:getSeqHelper', {}); });
}, {tier: 'api', description: 'awaits Bio\'s init through the platform (a call of a Bio service getter returns once initBio has run); not in the video'});

/** What the column's Renderer Settings > Font Size does, off camera: the monomers of a sequence column drawn larger. */
export const monomersDrawnAt = Given('the monomers of the {string} column are drawn at {int} px',
  async (page: Page, column: string, size: number) => {
    silent(page);
    await page.evaluate(([name, px]) => {
      const col = grok.shell.t?.col(name);
      if (!col)
        throw new Error(`the current table has no "${name}" column`);
      col.temp['.mm.cellRenderer.fontSize'] = px;
      col.temp['.mm.cellRenderer.settingsChanged'] = '1';
      col.dataFrame.fireValuesChanged();
    }, [column, size] as const);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'the renderer font size of a sequence column, as its Renderer Settings set it'});

/** The current row is drawn green, which tints what a clip colors by monomer. */
export const noRowCurrent = Given('no row of the table is current', async (page: Page) => {
  silent(page);
  await page.evaluate(() => { grok.shell.tv.dataFrame.currentRowIdx = -1; });
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx), {message: 'the current row'}).toBe(-1);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a clip leaves out: the current row of the current view\'s table cleared'});
