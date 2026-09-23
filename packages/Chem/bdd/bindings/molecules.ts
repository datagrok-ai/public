/* What the molecules in a column are, read through RDKit in the page (Chem:getRdKitModule): the
   same molecule in two notations, which rows a transform changed, fragments and substructures.
   Empty cells are compared as empty; a cell RDKit cannot read fails the step with its row. */
import {Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, gestures, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

type Pair = {row: number; a: string; b: string};

/** Canonical SMILES of every row of the named columns of the current table, '' for an empty cell;
 * flat drops the stereo marks (chiral tags, double-bond directions) before canonicalizing. */
async function canonical(page: Page, columns: string[], flat = false): Promise<{values: string[][]; unreadable: string[]}> {
  return page.evaluate(async ([names, noStereo]) => {
    const df = grok.shell.t;
    for (const n of names)
      if (!df.col(n))
        throw new Error(`no "${n}" column in ${df.name}; it has: ${df.columns.names().join(', ')}`);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const unreadable: string[] = [];
    const values = names.map((n: string) => {
      const col = df.col(n);
      const out: string[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const v = col.isNone(i) ? '' : String(col.get(i));
        if (v === '') {
          out.push('');
          continue;
        }
        let mol = null;
        try {
          mol = rdkit.get_mol(v);
          const smiles = mol.get_smiles();
          if (!noStereo)
            out.push(smiles);
          else {
            const bare = rdkit.get_mol(smiles.replace(/[@\/\\]/g, ''));
            out.push(bare.get_smiles());
            bare.delete();
          }
        }
        catch {
          out.push('');
          unreadable.push(`${n} row ${i + 1}`);
        }
        finally {
          mol?.delete();
        }
      }
      return out;
    });
    return {values, unreadable};
  }, [columns, flat] as [string[], boolean]);
}

function pairs(a: string[], b: string[]): Pair[] {
  return a.map((x, i) => ({row: i + 1, a: x, b: b[i]}));
}

async function expectSameMolecules(page: Page, x: string, y: string, flat: boolean): Promise<void> {
  const {values: [a, b], unreadable} = await canonical(page, [x, y], flat);
  expect(unreadable, 'cells RDKit could not read').toEqual([]);
  expect(a.filter(Boolean).length, `molecules in "${x}"`).toBeGreaterThan(0);
  const bad = pairs(a, b).filter((p) => p.a !== p.b).map((p) => `row ${p.row}: ${p.a} vs ${p.b}`);
  expect(bad.slice(0, 10), `rows where "${x}" and "${y}" hold different molecules (${bad.length})`).toEqual([]);
}

export const sameMolecules = Then('every molecule of {string} column should be the same as in {string} column', (page: Page, x: string, y: string) =>
  expectSameMolecules(page, x, y, false), {description: 'canonical SMILES of the two columns equal row by row: the same molecules in any notation'});

export const sameFlatMolecules = Then('every molecule of {string} column should be the same as in {string} column ignoring stereochemistry', (page: Page, x: string, y: string) =>
  expectSameMolecules(page, x, y, true), {description: 'canonical SMILES without chiral tags and double-bond directions equal row by row: the same atoms and bonds'});

export const noSameMolecule = Then('no molecule of {string} column should be the same as in {string} column', async (page: Page, x: string, y: string) => {
  const {values: [a, b], unreadable} = await canonical(page, [x, y]);
  expect(unreadable, 'cells RDKit could not read').toEqual([]);
  const same = pairs(a, b).filter((p) => p.a !== '' && p.a === p.b).map((p) => `row ${p.row}: ${p.a}`);
  expect(same, `rows where "${x}" repeats the molecule of "${y}"`).toEqual([]);
}, {description: 'no filled cell holds the molecule of the other column\'s cell in its row'});

export const changedMolecules = Then('{int} molecules of {string} column should differ from {string} column', async (page: Page, count: number, x: string, y: string) => {
  const {values: [a, b], unreadable} = await canonical(page, [x, y]);
  expect(unreadable, 'cells RDKit could not read').toEqual([]);
  const changed = pairs(a, b).filter((p) => p.a !== p.b).map((p) => p.row);
  expect(changed.length, `rows where "${x}" holds another molecule than "${y}": ${changed.join(', ')}`).toBe(count);
}, {description: 'the number of rows whose canonical SMILES differ between the two columns'});

export const someMoleculesChanged = Then('some but not all molecules of {string} column should differ from {string} column', async (page: Page, x: string, y: string) => {
  const {values: [a, b], unreadable} = await canonical(page, [x, y]);
  expect(unreadable, 'cells RDKit could not read').toEqual([]);
  const changed = pairs(a, b).filter((p) => p.a !== p.b).length;
  expect(changed, `rows where "${x}" holds another molecule than "${y}"`).toBeGreaterThan(0);
  expect(changed, `rows where "${x}" holds another molecule than "${y}"`).toBeLessThan(a.length);
});

/** The filled rows (1-based) whose molecule contains the query, an RDKit substructure match. */
function rowsContaining(page: Page, column: string, query: string): Promise<string[]> {
  return page.evaluate(async ([n, q]) => {
    const col = grok.shell.t.col(n);
    if (!col)
      throw new Error(`no "${n}" column in ${grok.shell.t.name}`);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const qmol = rdkit.get_qmol(q);
    const out: string[] = [];
    try {
      for (let i = 0; i < col.length; i++) {
        if (col.isNone(i) || col.get(i) === '')
          continue;
        const mol = rdkit.get_mol(col.get(i));
        try {
          if (mol.get_substruct_match(qmol) !== '{}')
            out.push(`row ${i + 1}`);
        }
        finally {
          mol.delete();
        }
      }
    }
    finally {
      qmol.delete();
    }
    return out;
  }, [column, query] as [string, string]);
}

export const noSubstructure = Then('no molecule of {string} column should contain {string}', async (page: Page, column: string, query: string) => {
  expect(await rowsContaining(page, column, query), `molecules of "${column}" containing ${query}`).toEqual([]);
}, {description: 'an RDKit substructure query (SMILES or SMARTS) matched against every filled cell'});

export const moleculesMatching = Then('{int} molecules of {string} column should contain {string}', async (page: Page, count: number, column: string, query: string) => {
  expect((await rowsContaining(page, column, query)).length, `molecules of "${column}" containing ${query}`).toBe(count);
});

export const noMoreFragments = Then('no molecule of {string} column should have more fragments than in {string} column', async (page: Page, x: string, y: string) => {
  const {values: [a, b]} = await canonical(page, [x, y]);
  const parts = (s: string) => s === '' ? 0 : s.split('.').length;
  const bad = pairs(a, b).filter((p) => parts(p.a) > parts(p.b)).map((p) => `row ${p.row}: ${p.a}`);
  expect(bad, `rows where "${x}" has more fragments than "${y}"`).toEqual([]);
}, {description: 'dot-separated components of the canonical SMILES, row by row'});

export const clipboardMolecule = Then('the clipboard should hold the molecule of row {int} of {string} column', async (page: Page, row: number, column: string) => {
  const text = await gestures.readClipboard(page);
  const r = await page.evaluate(async ([i, n, text]) => {
    const col = grok.shell.t.col(n);
    if (!col)
      throw new Error(`no "${n}" column in ${grok.shell.t.name}`);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const smiles = (v: string) => {
      const mol = rdkit.get_mol(v);
      try {
        return mol.get_smiles();
      }
      finally {
        mol.delete();
      }
    };
    let copied = '';
    try {
      copied = smiles(text);
    }
    catch {
      copied = `unreadable: ${text.slice(0, 60)}`;
    }
    return {copied, cell: smiles(String(col.get(i - 1)))};
  }, [row, column, text] as [number, string, string]);
  expect(r.copied, `the clipboard's molecule against row ${row} of "${column}"`).toBe(r.cell);
}, {description: 'the clipboard text read by RDKit (SMILES, SMARTS or molfile) is the cell\'s molecule'});

export const sortedBySimilarity = Then('the first {int} rows of grid should be in falling similarity to row {int} of {string} column', async (page: Page, count: number, row: number, column: string) => {
  await expect.poll(() => page.evaluate(() => (grok.shell.tv.grid.getRowOrder() as Int32Array)[0] + 1),
    {message: 'the table row the grid shows first'}).toBe(row);
  const r = await page.evaluate(async ([k, q, n]) => {
    const grid = grok.shell.tv.grid;
    const col = grid.dataFrame.col(n);
    if (!col)
      throw new Error(`no "${n}" column in ${grid.dataFrame.name}`);
    const order: number[] = Array.from(grid.getRowOrder() as Int32Array).slice(0, k);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const fp = (v: string) => {
      const mol = rdkit.get_mol(v);
      try {
        return mol.get_morgan_fp_as_uint8array(JSON.stringify({radius: 2, nBits: 2048}));
      }
      finally {
        mol.delete();
      }
    };
    const bits = (a: Uint8Array) => a.reduce((s, x) => s + x.toString(2).split('1').length - 1, 0);
    const query = fp(String(col.get(q - 1)));
    const scores = order.map((i) => {
      const f = fp(String(col.get(i)));
      const both = bits(f.map((x: number, j: number) => x & query[j]));
      const any = bits(f.map((x: number, j: number) => x | query[j]));
      return any === 0 ? 0 : both / any;
    });
    return {first: order[0] + 1, scores};
  }, [count, row, column] as [number, number, string]);
  expect(r.first, 'the table row the grid shows first').toBe(row);
  expect(r.scores[0], 'the similarity of the first row').toBe(1);
  for (let i = 1; i < r.scores.length; i++)
    expect(r.scores[i], `the similarity of grid row ${i + 1} against grid row ${i} (${r.scores.join(', ')})`).toBeLessThanOrEqual(r.scores[i - 1] + 1e-9);
}, {description: 'the grid\'s row order (getRowOrder): the query row first, then Tanimoto on Morgan fingerprints (radius 2, 2048 bits) never rising'});

/** The rows of the current table where the filter and an RDKit substructure match of the query
 * disagree, with how many rows contain the query, how many pass, and how many RDKit could not read. */
function filterAgainstQuery(page: Page, column: string, query: string):
  Promise<{wrong: string[]; matched: number; unreadable: number; passed: number}> {
  return page.evaluate(async ([c, q]) => {
    const df = grok.shell.t;
    const col = df.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${df.name}`);
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const qmol = rdkit.get_qmol(q);
    const wrong: string[] = [];
    let matched = 0;
    let unreadable = 0;
    try {
      for (let i = 0; i < df.rowCount; i++) {
        let hit = false;
        try {
          const mol = rdkit.get_mol(col.get(i));
          hit = mol.get_substruct_match(qmol) !== '{}';
          mol.delete();
        }
        catch {
          unreadable++;
        }
        if (hit)
          matched++;
        if (hit !== df.filter.get(i))
          wrong.push(`row ${i + 1} ${hit ? 'contains it but is filtered out' : 'passes without it'}`);
      }
    }
    finally {
      qmol.delete();
    }
    return {wrong, matched, unreadable, passed: df.filter.trueCount};
  }, [column, query] as [string, string]);
}

export const filterPassesMatching = Then('the filter should pass exactly the molecules of {string} column containing {string}', async (page: Page, column: string, query: string) => {
  const r = await filterAgainstQuery(page, column, query);
  expect(r.unreadable, `molecules of "${column}" RDKit could not read`).toBe(0);
  expect(r.wrong.slice(0, 10), `rows the filter disagrees with ${query} on (${r.wrong.length}; ${r.matched} contain it, ${r.passed} pass)`).toEqual([]);
}, {description: 'every row passes the filter if and only if its molecule contains the query (RDKit substructure match); polls nothing, so it follows a row-count step'});

export const readingIsRowMolecule = Then('the {string} reading of filter panel should be the molecule of row {int} of {string} column', async (page: Page, reading: string, row: number, column: string) => {
  await expect.poll(() => page.evaluate(async ([name, i, c]) => {
    // a cloned view leaves the first view's panel in the DOM, hidden
    const host = Array.from(document.querySelectorAll('[name="viewer-Filters"]')).find((e) => (e as HTMLElement).offsetParent !== null);
    const values = host == null ? {} : (window as any).__bdd.viewerOf(host).getWidgetStatus()?.values ?? {};
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const smiles = (v: string) => {
      if (!v)
        return '';
      const mol = rdkit.get_mol(v);
      try {
        return mol.get_smiles();
      }
      finally {
        mol.delete();
      }
    };
    return smiles(String(values[name] ?? '')) === smiles(String(grok.shell.t.col(c).get(i - 1)));
  }, [reading, row, column] as [string, number, string]), {message: `the "${reading}" reading as the molecule of row ${row} of "${column}"`}).toBe(true);
}, {description: 'the reading and the cell read by RDKit as the same molecule'});

export const cardsFromFilteredRows = Then('every card of {widget} should show a row that passes the filter', async (page: Page, target: ElementRef) => {
  let seen = '';
  await expect.poll(async () => {
    const r = await viewers.onViewer(page, target, (e) => {
      const v = (window as any).__bdd.viewerOf(e);
      const rows = String(v.getWidgetStatus()?.values?.['card row set'] ?? '').split(', ').filter(Boolean).map(Number);
      const df = grok.shell.t;
      return {rows, out: rows.filter((i: number) => !df.filter.get(i)).length, passing: df.filter.trueCount, limit: v.limit};
    });
    seen = `cards ${r.rows.join(', ')}; ${r.out} of them filtered out; ${r.passing} rows pass, limit ${r.limit}`;
    return r.out === 0 && r.rows.length === Math.min(r.passing, r.limit);
  }, {message: 'the cards against the rows that pass the filter'}).toBe(true).catch(() => {
    throw new Error(`the cards are not the rows that pass the filter: ${seen}`);
  });
}, {description: 'the viewer\'s "card row set" reading against the table filter: each card is a passing row, and there are as many cards as passing rows up to the limit'});

export const matrixDiagonal = Then('the similarity columns of table {string} should be symmetric, read {float} on the diagonal and less somewhere off it', async (page: Page, name: string, self: number) => {
  const r = await page.evaluate(([n, s]) => {
    const df = (grok.shell.tables ?? []).find((t: any) => t.name === n);
    if (!df)
      throw new Error(`no "${n}" table among ${(grok.shell.tables ?? []).map((t: any) => t.name).join(', ')}`);
    const cols = df.columns.toList().filter((c: any) => c.type === 'double' || c.type === 'float');
    const diagonal: number[] = [];
    let offDiagonal = 0;
    let asymmetric = 0;
    for (let i = 0; i < df.rowCount && i < cols.length; i++) {
      diagonal.push(cols[i].get(i));
      for (let j = 0; j < cols.length; j++) {
        if (j === i)
          continue;
        if (cols[j].get(i) < s)
          offDiagonal++;
        if (Math.abs(cols[j].get(i) - cols[i].get(j)) > 1e-5)
          asymmetric++;
      }
    }
    return {rows: df.rowCount, cols: cols.length, wrong: diagonal.filter((v) => Math.abs(v - s) > 1e-6).length, offDiagonal, asymmetric};
  }, [name, self] as [string, number]);
  expect(r.cols, `the similarity columns of "${name}" against its ${r.rows} rows`).toBe(r.rows);
  expect(r.wrong, `diagonal cells of "${name}" that are not ${self}`).toBe(0);
  expect(r.offDiagonal, `off-diagonal cells of "${name}" below ${self}`).toBeGreaterThan(0);
  expect(r.asymmetric, `pairs of "${name}" whose two cells differ`).toBe(0);
}, {description: 'a pairwise similarity table: one numeric column per row, every molecule fully similar to itself, and some pair less so'});

export const readingIsMolecule = Then('the {string} reading of {widget} should be the molecule {string}', async (page: Page, reading: string, target: ElementRef, smiles: string) => {
  let seen = '';
  await expect.poll(async () => {
    const r = await viewers.onViewer(page, target, async (e, arg: any) => {
      const value = String((window as any).__bdd.viewerOf(e).getWidgetStatus()?.values?.[arg.reading] ?? '');
      const rdkit = await grok.functions.call('Chem:getRdKitModule');
      const canonical = (v: string) => {
        if (!v)
          return '';
        const mol = rdkit.get_mol(v);
        try {
          return mol.get_smiles();
        }
        finally {
          mol.delete();
        }
      };
      return {value: canonical(value), want: canonical(arg.smiles)};
    }, {reading, smiles});
    seen = r.value;
    return r.value === r.want;
  }, {message: `the "${reading}" reading as the molecule ${smiles}`}).toBe(true).catch(() => {
    throw new Error(`the "${reading}" reading holds ${seen || 'nothing'}, not ${smiles}`);
  });
}, {description: 'a reading that holds a molecule in any notation, read by RDKit and compared with the one named'});

export const filterMatchesReading = Then('the filter should pass exactly the molecules of {string} column containing the {string} reading of {widget}', async (page: Page, column: string, reading: string, target: ElementRef) => {
  const query = await viewers.onViewer(page, target, (e, name: any) =>
    String((window as any).__bdd.viewerOf(e).getWidgetStatus()?.values?.[name] ?? ''), reading);
  expect(query, `the "${reading}" reading to filter by`).not.toBe('');
  const r = await filterAgainstQuery(page, column, query);
  expect(r.unreadable, `molecules of "${column}" RDKit could not read`).toBe(0);
  expect(r.matched, 'molecules containing the scaffold the reading holds').toBeGreaterThan(0);
  expect(r.wrong.slice(0, 10), `rows the filter disagrees with the "${reading}" reading on (${r.wrong.length}; ${r.matched} contain it, ${r.passed} pass)`).toEqual([]);
}, {description: 'the scaffold a viewer reports, matched against the column by RDKit, against what the table filter keeps'});
