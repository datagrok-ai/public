/* ---
sub_features_covered: [hmm.tour.screenshots]
--- */
/** A guided tour of the HMM package with Datagrok's own tutorial hints (a
 * popup that explains, a pulsing marker on the element) and a screenshot per
 * step, written to playwright/test-output/tour/. Not part of the regular
 * suite; run with:
 *   HMM_TOUR=1 grok test --host localhost --skip-build --skip-puppeteer --category tour */
import {test, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use({...specTestOptions, viewport: {width: 1440, height: 900}});

const OUT = 'playwright/test-output/tour';
const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';

async function openBioMenu(page: Page, group: string, leaf: string): Promise<void> {
  const groupLoc = page.locator(`[name="div-Bio---${group}"]`);
  if (!await groupLoc.isVisible().catch(() => false))
    await page.locator('[name="div-Bio"]').click();
  await groupLoc.waitFor({state: 'visible', timeout: 15_000});
  const box = (await groupLoc.boundingBox())!;
  await page.mouse.move(box.x + 5, box.y + box.height / 2, {steps: 12});
  await page.mouse.move(box.x + box.width - 15, box.y + box.height / 2, {steps: 12});
  const leafLoc = page.locator(`[name="div-Bio---${group}---${leaf}"]`);
  await leafLoc.waitFor({state: 'visible', timeout: 15_000});
  await leafLoc.click();
  await page.mouse.move(5, 450);
}

async function setChoice(page: Page, dialog: string, input: string, value: string): Promise<void> {
  await page.locator(`${dialog} [name="input-host-${input}"] select`).selectOption({label: value});
}

/** Datagrok's tutorial hint next to an element (CSS selector) or over a grid
 * column header ({column}): an explaining popup and a pulsing marker. */
async function explain(page: Page, target: string | {column: string}, title: string, text: string,
  position: 'top' | 'bottom' | 'left' | 'right' = 'right'): Promise<void> {
  await page.evaluate(({target, title, text, position}) => {
    let el: HTMLElement;
    if (typeof target === 'string')
      el = document.querySelector(target) as HTMLElement;
    else {
      const grid = grok.shell.tv.grid;
      const col = grid.columns.byName(target.column)!;
      const r = grid.root.getBoundingClientRect();
      el = document.createElement('div');
      el.className = 'tour-anchor';
      Object.assign(el.style, {position: 'fixed', left: `${r.left + col.left}px`, top: `${r.top}px`,
        width: `${Math.min(col.width, r.width - col.left)}px`, height: `${grid.colHeaderHeight}px`,
        pointerEvents: 'none', zIndex: '5'});
      document.body.appendChild(el);
    }
    const content = ui.divV([ui.h3(title), ui.divText(text)], {style: {maxWidth: '320px', lineHeight: '1.45'}});
    const hint = ui.hints.addHint(el, content, position);
    hint.classList.add('ui-hint-popup-green-1');
    ui.hints.addHintIndicator(el, false);
  }, {target, title, text, position});
  await page.waitForTimeout(800);
}

async function clearHints(page: Page): Promise<void> {
  await page.evaluate(() => {
    document.querySelectorAll('.ui-hint-popup, .ui-hint-blob, .ui-hint-overlay, .tour-anchor')
      .forEach((e) => e.remove());
    document.querySelectorAll('.ui-hint-target').forEach((e) => e.classList.remove('ui-hint-target'));
  });
}

/** A screenshot of the whole window, or of the box around `selectors` (and the hints). */
async function shot(page: Page, name: string, selectors: string[] = []): Promise<void> {
  // Toasts and the local server's "debugging packages" balloon are not part of the feature.
  await page.evaluate(() => document.querySelectorAll('.d4-balloon').forEach((e) => e.remove()));
  if (selectors.length === 0) {
    await page.screenshot({path: `${OUT}/${name}.png`});
    return;
  }
  const box = await page.evaluate((sels) => {
    const rects = [...sels, '.ui-hint-popup'].flatMap((s) => Array.from(document.querySelectorAll(s)))
      .map((e) => e.getBoundingClientRect()).filter((r) => r.width > 0);
    const x = Math.max(0, Math.min(...rects.map((r) => r.left)) - 24);
    const y = Math.max(0, Math.min(...rects.map((r) => r.top)) - 24);
    const right = Math.min(window.innerWidth, Math.max(...rects.map((r) => r.right)) + 24);
    const bottom = Math.min(window.innerHeight, Math.max(...rects.map((r) => r.bottom)) + 24);
    return {x, y, width: right - x, height: bottom - y};
  }, selectors);
  await page.screenshot({path: `${OUT}/${name}.png`, clip: box});
}

/** Put `names` first in the grid (in that order), widened as given. */
async function showColumns(page: Page, columns: [string, number?][]): Promise<void> {
  await page.evaluate((cols) => {
    const grid = grok.shell.tv.grid;
    const names = cols.map(([n]) => n);
    const rest = grok.shell.tv.dataFrame.columns.names().filter((n) => !names.includes(n) && !n.startsWith('~'));
    grid.columns.setOrder([...names, ...rest]);
    for (const [n, w] of cols) if (w) grid.columns.byName(n)!.width = w;
    grid.scrollToCell(names[0], 0);
  }, columns);
  await page.waitForTimeout(800);
}

test('HMM package tour: explained screenshots', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async ({path}) => {
    document.body.classList.add('selenium');
    grok.shell.closeAll();
    const src = await grok.dapi.files.readCsv(path);
    // The first 40 chains without a signal peptide (leaders would only widen the aligned views).
    const hc = src.col('AntibodyHC')!;
    const keep = Array.from({length: src.rowCount}, (_, i) => i).filter((i) => !(hc.get(i) ?? '').startsWith('M'))
      .slice(0, 40);
    const df = src.clone(DG.BitSet.create(src.rowCount, (i) => keep.includes(i)));
    df.name = 'antibodies';
    grok.shell.addTableView(df);
    for (let i = 0; i < 200; i++) {
      if (df.col('AntibodyHC')?.semType === 'Macromolecule' && document.querySelector('[name="viewer-Grid"] canvas'))
        break;
      await new Promise((r) => setTimeout(r, 200));
    }
    grok.shell.windows.simpleMode = true;
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showBrowse = false;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showHelp = false;
  }, {path: ANTIBODIES});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 60_000});

  const numbering = '[name="dialog-Apply-Antibody-Numbering"]';
  await softStep('01 numbering dialog with engine options', async () => {
    await openBioMenu(page, 'Annotate', 'Apply-Numbering-Scheme...');
    await page.locator(numbering).waitFor({timeout: 60_000});
    await setChoice(page, numbering, 'Engine', 'ANARCI (HMMER)');
    await setChoice(page, numbering, 'Scheme', 'imgt');
    await page.locator(`${numbering} [name="input-host-Domain"]`).waitFor({timeout: 30_000});
    await explain(page, `${numbering} [name="input-host-Species"]`, 'Engine options',
      'ANARCI (HMMER) adds its own options to Bio\'s numbering dialog: the preferred species of the germline ' +
      'models, the bit score a domain needs, and which domain of an scFv to number (first, heavy or light).');
    await shot(page, '01-numbering-dialog', [numbering]);
    await clearHints(page);
    await page.locator(`${numbering} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC (aligned)'), null,
      {timeout: 120_000});
  });

  await softStep('02 numbered and aligned', async () => {
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 360], ['AntibodyHC (aligned)', 720]]);
    await explain(page, {column: 'AntibodyHC (aligned)'}, 'IMGT numbering in the browser',
      'Every chain is numbered by ANARCI\'s method on the HMMER engine and aligned by IMGT position. ' +
      'Framework and CDR regions are drawn on both columns. Identical to unmodified ANARCI.', 'bottom');
    await shot(page, '02-numbering-result');
    await clearHints(page);
  });

  await softStep('03 germlines and species', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      await grok.functions.call('HMM:antibodyGermlines', {table: df, sequence: df.col('AntibodyHC'),
        species: 'human and mouse'});
    });
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 300], ['AntibodyHC chain', 80],
      ['AntibodyHC species', 80], ['AntibodyHC v gene', 120], ['AntibodyHC v identity', 90],
      ['AntibodyHC j gene', 100], ['AntibodyHC evalue', 100]]);
    await explain(page, {column: 'AntibodyHC v gene'}, 'Germlines and species',
      'Bio › Annotate › Germlines and Species (ANARCI) adds chain type, species, the closest V and J genes ' +
      'with their identities, and the HMMER E-value of each chain.', 'bottom');
    await shot(page, '03-germlines');
    await clearHints(page);
  });


  await softStep('04 context panels: ANARCI and Pfam domains', async () => {
    await page.evaluate(() => {
      grok.shell.windows.showContextPanel = true;
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, 'AntibodyHC'));
    });
    await page.evaluate(async () => {
      const header = async (text: string) => {
        for (let t = 0; t < 150; t++) {
          const h = Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
            .find((e) => (e.textContent ?? '').trim() === text) as HTMLElement | undefined;
          if (h) return h;
          await new Promise((r) => setTimeout(r, 200));
        }
        return null;
      };
      const expand = async (text: string) => {
        const h = await header(text);
        if (!h) return;
        const content = h.nextElementSibling as HTMLElement | null;
        if (!content || content.offsetHeight === 0) h.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
      };
      await expand('Bioinformatics');
      await expand('ANARCI');
      await expand('Pfam Domains');
      for (let t = 0; t < 100; t++) {
        const text = document.querySelector('.grok-prop-panel')?.textContent ?? '';
        if (text.includes('V gene') && text.includes('C1-set')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    });
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 420]]);
    const anchor = await page.evaluate(() => {
      const h = Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
        .find((e) => (e.textContent ?? '').trim() === 'ANARCI') as HTMLElement | undefined;
      if (h) h.setAttribute('data-tour', 'anarci-panel');
      return !!h;
    });
    if (anchor) {
      await explain(page, '[data-tour="anarci-panel"]', 'Panels for one sequence',
        'Click a sequence: the context panel shows its ANARCI result (chain, species, germlines, E-value, IMGT ' +
        'CDRs) and its Pfam domains, computed in a few milliseconds.', 'left');
    }
    await shot(page, '04-context-panels');
    await clearHints(page);
    await page.evaluate(() => { grok.shell.windows.showContextPanel = false; });
  });

  const domains = '[name^="dialog-Find-Domains"]';
  await softStep('04 find domains dialog', async () => {
    await openBioMenu(page, 'Annotate', 'Find-Domains-(HMMER)...');
    await page.locator(domains).waitFor({timeout: 60_000});
    await explain(page, `${domains} [name="input-host-Clan-competition"]`, 'Clan competition',
      'Overlapping matches of related families (one Pfam clan, such as the immunoglobulin-like domains) are ' +
      'resolved as Pfam does: the most significant one is kept. Identical to Pfam\'s pfam_scan.pl.');
    await shot(page, '05-domains-dialog', [domains]);
    // The note that HMMER is reimplemented here is in the tooltips, not in the dialog.
    await clearHints(page);
    // Hovering fills the input's tooltip; automated browsers keep tooltips hidden, so the
    // text Datagrok put in it is shown in a tour hint instead.
    const label = page.locator(`${domains} [name="input-host-Library"] label`).first();
    const box = (await label.boundingBox())!;
    await page.mouse.move(box.x + 4, box.y + box.height / 2);
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2, {steps: 5});
    const tooltip = await page.waitForFunction(() => Array.from(document.querySelectorAll('.d4-tooltip'))
      .map((e) => e.textContent ?? '').find((t) => /reimplemented/.test(t)), null, {timeout: 8_000})
      .then((h) => h.jsonValue() as Promise<string>, () => null);
    await page.mouse.move(5, 5);
    if (tooltip) {
      await explain(page, `${domains} [name="input-host-Library"]`, 'Tooltip of Library', tooltip);
      await shot(page, '05b-tooltip-note', [domains]);
      await clearHints(page);
    }
    await clearHints(page);
    await page.locator(`${domains} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC domains'), null,
      {timeout: 120_000});
  });

  await softStep('05 domains on the sequences', async () => {
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 520], ['AntibodyHC domains', 220]]);
    await explain(page, {column: 'AntibodyHC domains'}, 'Pfam domains',
      'Each Fab heavy chain shows one V-set and one C1-set domain from the built-in Pfam library (61 families ' +
      'for biologics), with Pfam\'s gathering thresholds. Domains are colored on the sequences too.', 'bottom');
    await shot(page, '06-domains-result');
    await clearHints(page);
  });


  await softStep('07 clan competition on the light chains', async () => {
    const counts = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const lc = df.col('AntibodyLC')!;
      const summary = async (clanCompetition: boolean) => {
        const found: DG.DataFrame = await grok.functions.call('HMM:findDomains', {table: df, sequence: lc,
          library: 'builtin', evalue: 0.01, gathering: true, clanCompetition});
        const text: [number, string][][] = Array.from({length: df.rowCount}, () => []);
        for (let i = 0; i < found.rowCount; i++) {
          text[found.get('row', i)].push([found.get('from', i),
            `${found.get('model', i)} ${found.get('from', i)}-${found.get('to', i)}`]);
        }
        return text.map((t) => t.sort((a, b) => a[0] - b[0]).map((x) => x[1]).join('; '));
      };
      const all = await summary(false);
      const kept = await summary(true);
      df.columns.addNewString('AntibodyLC: all Pfam matches').init((i) => all[i]);
      df.columns.addNewString('AntibodyLC: clan competition').init((i) => kept[i]);
      const count = (texts: string[]) => texts.reduce((n, t) => n + (t ? t.split('; ').length : 0), 0);
      return {all: count(all), kept: count(kept)};
    });
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyLC', 300], ['AntibodyLC: all Pfam matches', 420],
      ['AntibodyLC: clan competition', 260]]);
    await explain(page, {column: 'AntibodyLC: clan competition'}, 'Clan competition at work',
      'In light chains, I-set, Ig_3 and ig all match the same stretch as V-set. They belong to one Pfam clan, ' +
      `so only the most significant match is kept, as Pfam reports it (${counts.all - counts.kept} of ${counts.all} ` +
      'matches dropped here).', 'bottom');
    await shot(page, '07-clan-competition');
    await clearHints(page);
  });

  const profile = '[name="dialog-Profile-HMM-Search"]';
  await softStep('08 profile HMM search dialog', async () => {
    await openBioMenu(page, 'Search', 'Profile-HMM-Search...');
    await page.locator(profile).waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const select = document.querySelector('[name="dialog-Profile-HMM-Search"] [name="input-host-Family"] select');
      return !!select && (select as HTMLSelectElement).options.length > 1;
    }, null, {timeout: 60_000});
    await explain(page, `${profile} [name="input-host-Family"] select`, 'Search the table with one family',
      'hmmsearch asks every row whether it belongs to a family: a built-in Pfam family, any Pfam accession ' +
      '(fetched from InterPro and cached in the browser) or an HMM file. The table is the database.', 'bottom');
    await shot(page, '08-profile-search-dialog', [profile]);
    await clearHints(page);
    await page.locator(`${profile} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('V-set E-value'), null,
      {timeout: 120_000});
  });

  await softStep('09 profile HMM search result', async () => {
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 420], ['V-set score', 100], ['V-set E-value', 110],
      ['V-set domains', 80], ['V-set significant', 90], ['V-set regions', 100]]);
    await explain(page, {column: 'V-set E-value'}, 'Scores for every row',
      'Bit score, E-value (the column is the database), domain count, significance and the matched regions, ' +
      'identical to C hmmsearch even though the search is split across workers.', 'bottom');
    await shot(page, '09-profile-search-result');
    await clearHints(page);
  });

  const msa = '[name="dialog-MSA"]';
  await softStep('06 align to profile in the MSA dialog', async () => {
    await openBioMenu(page, 'Analyze', 'MSA...');
    await page.locator(msa).waitFor({timeout: 60_000});
    await setChoice(page, msa, 'Engine', 'Profile HMM (hmmalign)');
    await page.locator(`${msa} [name="input-host-Profile"] select`).waitFor({timeout: 30_000});
    await setChoice(page, msa, 'Profile', 'Antibody/TCR germline (auto)');
    await page.locator(`${msa} [name="input-host-Trim"] input[type="checkbox"]`).check();
    await explain(page, `${msa} [name="input-host-Engine"] select`, 'Align to a profile',
      'Protein and nucleotide columns can now use registered engines next to Kalign. Profile HMM (hmmalign) ' +
      'aligns every row to one model: a Pfam family, the best ANARCI germline model or an HMM file.');
    await shot(page, '10-msa-dialog', [msa]);
    await clearHints(page);
    await page.locator(`${msa} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.names().some((n) => n.startsWith('msa(')),
      null, {timeout: 120_000});
  });

  await softStep('07 profile alignment', async () => {
    const name = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names().find((n) => n.startsWith('msa('))!);
    await showColumns(page, [['Antibody_ID', 90], [name, 1000]]);
    await explain(page, {column: name}, 'Columns follow the model',
      'Each column is a position of the germline model (or an insertion between two positions), so alignments ' +
      'of different tables to one model line up; Trim drops residues outside the model. Identical to C ' +
      'hmmalign.', 'bottom');
    await shot(page, '11-msa-result');
    await clearHints(page);
  });

  const homology = '[name^="dialog-Homology-Search"]';
  await softStep('08 homology search dialog', async () => {
    await page.evaluate(() => { grok.shell.tv.dataFrame.currentRowIdx = 2; });
    await openBioMenu(page, 'Search', 'Homology-Search-(phmmer)...');
    await page.locator(homology).waitFor({timeout: 60_000});
    await explain(page, `${homology} [name="input-host-Query"]`, 'Homology search',
      'phmmer scores every row against one query (the current row or a typed sequence) with bit scores and ' +
      'E-values instead of percent identity. The table is the database.');
    await shot(page, '12-homology-dialog', [homology]);
    await clearHints(page);
    await page.locator(`${homology} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('phmmer row 3 E-value'), null,
      {timeout: 120_000});
  });

  await softStep('09 homology search result', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.rows.select((r) => false);
    });
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 420], ['phmmer row 3 score', 110],
      ['phmmer row 3 E-value', 120], ['phmmer row 3 significant', 90], ['phmmer row 3 regions', 120]]);
    await explain(page, {column: 'phmmer row 3 E-value'}, 'E-values for every row',
      'Row 3 itself scores best; other heavy chains follow by E-value. The query model, calibration and ' +
      'statistics are exactly phmmer\'s.', 'bottom');
    await shot(page, '13-homology-result');
    await clearHints(page);
  });

  const build = '[name="dialog-Build-Profile-HMM"]';
  await softStep('10 build a profile from selected rows', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i) => i < 15);
    });
    await openBioMenu(page, 'Search', 'Build-Profile-HMM...');
    await page.locator(build).waitFor({timeout: 60_000});
    await page.locator(`${build} [name="input-host-Name"] input`).fill('heavy15');
    await explain(page, `${build} [name="input-host-Rows"]`, 'Build a profile',
      'hmmbuild turns the 15 selected aligned rows into a profile HMM, searches the table with it, ' +
      'annotates the hits and opens the profile\'s logo in WebLogo. The model file is byte-identical to C ' +
      'hmmbuild\'s and is kept with the table.');
    await shot(page, '14-build-dialog', [build]);
    await clearHints(page);
    await page.locator(`${build} [name="button-OK"]`).click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('heavy15 E-value'), null,
      {timeout: 120_000});
    await page.waitForFunction(() => Array.from(grok.shell.tv.viewers).some((v) =>
      v.type === 'WebLogo' && (v as any).profile === 'heavy15' && (v as any).positions?.length > 0), null,
    {timeout: 60_000});
  });

  await softStep('11 search with the built profile, and its logo', async () => {
    await showColumns(page, [['Antibody_ID', 90], ['AntibodyHC', 260], ['heavy15 score', 100],
      ['heavy15 E-value', 110], ['heavy15 regions', 100]]);
    await explain(page, '[name="viewer-WebLogo"]', 'The profile\'s logo, in WebLogo',
      'Bio\'s WebLogo draws the profile: letter heights are hmmlogo\'s values (bits), so conserved framework ' +
      'positions stand tall and CDRs are low. Clicking a letter selects the rows of the alignment that have it.',
    'top');
    await shot(page, '15-build-result');
    await clearHints(page);
  });

  await softStep('12 logo of a Pfam family; the WebLogo Profile property', async () => {
    await page.evaluate(() => {
      for (const v of Array.from(grok.shell.tv.viewers)) if (v.type === 'WebLogo') v.close();
    });
    await openBioMenu(page, 'Analyze', 'Profile-HMM-Logo...');
    const logo = '[name="dialog-Profile-HMM-Logo"]';
    await page.locator(logo).waitFor({timeout: 60_000});
    await page.locator(`${logo} [name="input-host-Profile-from"] select`)
      .selectOption({label: 'Pfam family (built-in)'});
    await page.locator(`${logo} [name="input-host-Family"] select`).selectOption({label: 'C1-set'});
    await explain(page, `${logo} [name="input-host-Sequence"]`, 'Column and profile',
      'Choose the WebLogo\'s sequence column and a profile: a Pfam family of the built-in library, a Pfam ' +
      'accession, an ANARCI germline model or a profile built in this table. A column aligned to a profile ' +
      'proposes that profile.');
    await shot(page, '16-logo-dialog', [logo]);
    await clearHints(page);
    await page.locator(`${logo} [name="button-OK"]`).click();
    await page.waitForFunction(() => Array.from(grok.shell.tv.viewers).some((v) =>
      v.type === 'WebLogo' && (v as any).positions?.length === 86), null, {timeout: 60_000});
    // A taller strip, and the viewer's properties: Profile is an ordinary WebLogo setting.
    await page.evaluate(() => {
      const v = Array.from(grok.shell.tv.viewers).find((x) => x.type === 'WebLogo')!;
      grok.shell.tv.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'C1-set logo', 0.55);
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = v;
    });
    await page.waitForTimeout(1500);
    await explain(page, '[name="viewer-WebLogo"]', 'A Pfam family\'s logo',
      'C1-set from the built-in library. Change it in the Profile property (a family, an accession, ' +
      'germline:human_H or a profile built in the table); position width, Entropy or 100% heights and the ' +
      'range slider work as for any WebLogo.', 'top');
    await shot(page, '17-logo');
    await clearHints(page);
  });
  finishSpec(stepErrors);
});
