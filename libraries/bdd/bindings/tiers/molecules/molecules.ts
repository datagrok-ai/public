/* The `molecules` tier: what a molecule is, read through Chem's RDKit in the page (Chem:getRdKitModule), so a claim
   holds for the molecule in any notation: a reading of a widget that holds one (a sketcher's "smiles"), a cell of a
   column. And a table written in a feature made a table of molecules, as when a file opens. A feature using these needs
   Chem on the stand; a package other than Chem gates on it (`the "Chem" package is installed`). */
import {Page} from '@playwright/test';
import {expect} from '../../../src/runtime/patience.js';
import {Given, Then} from '../../../src/registry.js';
import type {ElementRef} from '../../../src/runtime/args.js';
import * as v from '../../../src/runtime/viewers.js';

declare const grok: any;

/** Canonical SMILES of the molecules, read by RDKit in the page; '' for an empty one, null for one RDKit cannot read. */
function canonicalOf(page: Page, molecules: string[]): Promise<(string | null)[]> {
  return page.evaluate(async (list) => {
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    return list.map((m: string) => {
      if (!m)
        return '';
      let mol = null;
      try {
        mol = rdkit.get_mol(m);
        return mol.get_smiles() as string;
      }
      catch {
        return null;
      }
      finally {
        mol?.delete();
      }
    });
  }, molecules);
}

/** A reading of the widget, as text ('' when it reports none). */
function readingOf(page: Page, target: ElementRef, reading: string): Promise<string> {
  return v.onViewer(page, target, (e, name: any) =>
    String((window as any).__bdd.viewerOf(e).getWidgetStatus()?.values?.[name] ?? ''), reading);
}

/** A cell of the current table as text, rows counted from 1. */
function cellOf(page: Page, column: string, row: number): Promise<string> {
  return page.evaluate(([c, r]) => {
    const df = grok.shell.t;
    const col = df?.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df?.name}; it has: ${df?.columns.names().join(', ')}`);
    return col.isNone(r - 1) ? '' : String(col.get(r - 1));
  }, [column, row] as [string, number]);
}

export const readingIsMolecule = Then('the {string} reading of {widget} should be the molecule {string}', async (page: Page, reading: string, target: ElementRef, smiles: string) => {
  let seen = '';
  await expect.poll(async () => {
    const value = await readingOf(page, target, reading);
    const [got, want] = await canonicalOf(page, [value, smiles]);
    seen = got ?? value;
    return got !== null && got === want;
  }, {message: `the "${reading}" reading as the molecule ${smiles}`}).toBe(true).catch(() => {
    throw new Error(`the "${reading}" reading holds ${seen || 'nothing'}, not ${smiles}`);
  });
}, {description: 'a reading that holds a molecule in any notation, read by RDKit and compared with the one named'});

export const readingIsRowMoleculeOf = Then('the {string} reading of {widget} should be the molecule in row {int} of {string} column',
  async (page: Page, reading: string, target: ElementRef, row: number, column: string) => {
    let seen = {value: '', cell: ''};
    await expect.poll(async () => {
      const [value, cell] = await canonicalOf(page, [await readingOf(page, target, reading), await cellOf(page, column, row)]);
      seen = {value: value ?? '', cell: cell ?? ''};
      return seen.value !== '' && seen.value === seen.cell;
    }, {message: `the "${reading}" reading as the molecule in row ${row} of "${column}"`}).toBe(true).catch(() => {
      throw new Error(`the "${reading}" reading holds ${seen.value || 'nothing'}, not ${seen.cell}`);
    });
  }, {description: 'a reading that holds a molecule, read by RDKit and compared with the cell\'s'});

export const rowMolecule = Then('the molecule in row {int} of {string} column should be {string}', async (page: Page, row: number, column: string, smiles: string) => {
  let seen = {cell: '', want: ''};
  await expect.poll(async () => {
    const [cell, want] = await canonicalOf(page, [await cellOf(page, column, row), smiles]);
    seen = {cell: cell ?? 'a value RDKit cannot read', want: want ?? ''};
    return cell !== null && cell === want;
  }, {message: `the molecule in row ${row} of "${column}"`}).toBe(true).catch(() => {
    throw new Error(`the molecule in row ${row} of "${column}" is ${seen.cell || 'empty'}, not ${seen.want}`);
  });
}, {description: 'the cell and the SMILES read by RDKit as the same molecule, stereochemistry included; polled, as a host writes the cell on its OK'});

export const detectTypes = Given('the semantic types of the current table are detected', async (page: Page) => {
  await page.evaluate(async () => {
    await grok.data.detectSemanticTypes(grok.shell.t);
  });
}, {tier: 'api', description: 'the platform\'s detectors over the current table, as when a file opens: a SMILES column becomes Molecule, units smiles'});
