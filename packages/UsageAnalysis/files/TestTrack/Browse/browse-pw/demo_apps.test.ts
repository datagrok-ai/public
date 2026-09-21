import { test, expect, Page } from '@playwright/test';
import { CONTEXT_PANEL, TREE_EXPAND_ARROW } from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  expandTreeGroup,
  watchErrors,
  expectNoErrors,
} from './helpers';

interface Demo {
  path: string[];
  view: string;
  tags?: Array<'heavy' | 'dashboard' | 'skip'>;
  ref?: string;

  opensForm?: boolean;

  envIgnore?: RegExp;
}

const DEMOS: Demo[] = [

  { path: ['Apps', 'Demo', 'Cheminformatics', 'Med Chem'], view: 'Med Chem' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Chemical Space'], view: 'Chemical Space', tags: ['skip'], ref: 'GROK-14320' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Molecule Activity Cliffs'], view: 'Molecule Activity Cliffs', tags: ['skip', 'dashboard'], ref: 'GROK-14320' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'R-Group Analysis'], view: 'R-Group Analysis', tags: ['skip', 'dashboard'], ref: 'GROK-14320' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Matched Molecular Pairs'], view: 'Matched Molecular Pairs' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Similarity & Diversity Search'], view: 'Similarity & Diversity Search' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Scaffold Tree'], view: 'Scaffold Tree' },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Database Queries'], view: 'Database Queries', tags: ['heavy'], opensForm: true, envIgnore: /Connection to .* refused|postmaster is accepting TCP\/IP connections|Stack trace \w+/i },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Admetica'], view: 'Admetica', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Cheminformatics', 'Retrosynthesis'], view: 'Retrosynthesis', tags: ['heavy'] },

  { path: ['Apps', 'Demo', 'Bioinformatics', 'Peptide SAR'], view: 'Peptide SAR', tags: ['dashboard', 'heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Sequence Activity Cliffs'], view: 'Sequence Activity Cliffs', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'siRNA'], view: 'siRNA', tags: ['skip', 'heavy'], ref: 'demoSkip' },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Antibodies'], view: 'Antibodies', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Sequence Space'], view: 'Sequence Space', tags: ['dashboard', 'heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Similarity, Diversity'], view: 'Similarity, Diversity', tags: ['heavy'], ref: 'GROK-18050' },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Atomic Level'], view: 'Atomic Level', tags: ['skip', 'heavy'], ref: 'demoSkip' },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Docking'], view: 'Docking', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Docking Conformations'], view: 'Docking Conformations', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Bioinformatics', 'Proteins'], view: 'Proteins', tags: ['heavy'] },

  { path: ['Apps', 'Demo', 'Data Access', 'Table Linking'], view: 'Table Linking' },
  { path: ['Apps', 'Demo', 'Data Access', 'Files'], view: 'Files' },
  { path: ['Apps', 'Demo', 'Data Access', 'Databases'], view: 'Databases', tags: ['heavy'] },

  { path: ['Apps', 'Demo', 'Visualization', 'Data Flow and Hierarchy', 'Network Diagram'], view: 'Network Diagram' },
  { path: ['Apps', 'Demo', 'Visualization', 'Data Flow and Hierarchy', 'Tree'], view: 'Tree' },
  { path: ['Apps', 'Demo', 'Visualization', 'Data Flow and Hierarchy', 'Tree Map'], view: 'Tree Map', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Visualization', 'Data Separation', 'Trellis Plot'], view: 'Trellis Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'Data Separation', 'Matrix Plot'], view: 'Matrix Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Scatter Plot'], view: 'Scatter Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Bar Chart'], view: 'Bar Chart' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Line Chart'], view: 'Line Chart' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Histogram'], view: 'Histogram' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Pie Chart'], view: 'Pie Chart' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', '3D Scatter Plot'], view: '3D Scatter Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Tile Viewer'], view: 'Tile Viewer' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Density Plot'], view: 'Density Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Filters'], view: 'Filters' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Heatmap'], view: 'Heatmap' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Markup'], view: 'Markup' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Radar'], view: 'Radar' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Sunburst'], view: 'Sunburst' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Chord'], view: 'Chord' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Sankey'], view: 'Sankey' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Surface Plot'], view: 'Surface Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Timelines'], view: 'Timelines' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Word Cloud'], view: 'Word Cloud' },
  { path: ['Apps', 'Demo', 'Visualization', 'General', 'Data Annotations'], view: 'Data Annotations' },
  { path: ['Apps', 'Demo', 'Visualization', 'Geographical', 'Map'], view: 'Map', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Visualization', 'Input and Edit', 'Grid'], view: 'Grid' },
  { path: ['Apps', 'Demo', 'Visualization', 'Input and Edit', 'Form'], view: 'Form' },
  { path: ['Apps', 'Demo', 'Visualization', 'Statistical', 'Box Plot'], view: 'Box Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'Statistical', 'Correlation Plot'], view: 'Correlation Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'Statistical', 'PC Plot'], view: 'PC Plot' },
  { path: ['Apps', 'Demo', 'Visualization', 'Statistical', 'Pivot Table'], view: 'Pivot Table' },
  { path: ['Apps', 'Demo', 'Visualization', 'Statistical', 'Statistics'], view: 'Statistics' },
  { path: ['Apps', 'Demo', 'Visualization', 'Time and Date', 'Calendar'], view: 'Calendar' },

  { path: ['Apps', 'Demo', 'Compute', 'Diff Studio'], view: 'Diff Studio', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Compute', 'Multivariate Analysis'], view: 'Multivariate Analysis', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Compute', 'PK-PD Modeling'], view: 'PK-PD Modeling', tags: ['heavy'] },
  { path: ['Apps', 'Demo', 'Compute', 'Bioreactor'], view: 'Bioreactor', tags: ['heavy'] },

  { path: ['Apps', 'Demo', 'Curves', 'Curve Fitting'], view: 'Curve Fitting' },
  { path: ['Apps', 'Demo', 'Curves', 'Assay Curves'], view: 'Assay Curves' },
];

function treeName(path: string[]): string {
  return 'tree-' + path.map((s) => s.replace(/ /g, '-')).join('---');
}

async function expandByTreeName(page: Page, name: string): Promise<void> {
  const node = page.locator(`[name="${name}"]`).first();
  await node.waitFor({ state: 'visible', timeout: 20_000 });
  await node.scrollIntoViewIfNeeded();
  const tri = node.locator(TREE_EXPAND_ARROW).first();
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false);
  if (expanded) return;
  if (await tri.isVisible().catch(() => false)) await tri.click();
  else await node.click();
  await page.waitForTimeout(900);
}

async function openDemo(page: Page, demo: Demo): Promise<void> {
  await ensureBrowsePanelOpen(page);
  await expandTreeGroup(page, 'Apps'); 

  for (let i = 2; i < demo.path.length; i++)
    await expandByTreeName(page, treeName(demo.path.slice(0, i)));

  const leaf = page.locator(`[name="${treeName(demo.path)}"]`).first();
  await leaf.waitFor({ state: 'visible', timeout: 20_000 });
  await leaf.scrollIntoViewIfNeeded();
  await leaf.click();
}

function norm(s: string): string {
  return s.toLowerCase().replace(/[^a-z0-9]+/g, '');
}

async function waitForDemoView(page: Page, demo: Demo, timeout: number): Promise<void> {
  const slug = norm(demo.path.slice(2).join('/'));
  const leaf = norm(demo.view);
  await expect
    .poll(async () => {
      const info = await page.evaluate(() => ({
        curPath: (window as any).grok.shell.v?.path ?? '',
        names: Array.from((window as any).grok.shell.views).map((v: any) => v.name ?? ''),
      }));
      if (norm(info.curPath).includes(slug)) return true;
      if (info.names.some((n: string) => {
        const a = norm(n);
        return !!a && (a.includes(leaf) || leaf.includes(a));
      })) return true;

      if (demo.opensForm) {
        const form = page.locator('.CodeMirror, .d4-dialog')
          .or(page.getByText(/Query SQL|Query Input Form/i)).first();
        if (await form.isVisible().catch(() => false)) return true;
      }
      return false;
    }, {
      message: `Demo view "${demo.view}" should open (path ~ ${demo.path.slice(2).join('/')})`,
      timeout,
      intervals: [500, 1000, 1500, 2000, 3000],
    })
    .toBe(true);
}

async function pokeView(page: Page): Promise<void> {
  const doc = page.locator('.d4-root .grok-view, .layout-workarea .d4-view, .document-manager .tab-content').first();
  if (await doc.isVisible().catch(() => false)) {
    const box = await doc.boundingBox().catch(() => null);
    if (box) await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2).catch(() => undefined);
  }

  const gridCanvas = page.locator('canvas.d4-grid, .d4-grid canvas').first();
  if (await gridCanvas.isVisible().catch(() => false)) {
    const gb = await gridCanvas.boundingBox().catch(() => null);
    if (gb) await page.mouse.click(gb.x + Math.min(60, gb.width / 2), gb.y + Math.min(80, gb.height / 2)).catch(() => undefined);
  }

  const viewer = page.locator('.d4-viewer:not(.d4-grid)').first();
  if (await viewer.isVisible().catch(() => false)) await viewer.hover().catch(() => undefined);
  await page.waitForTimeout(800);
}

test.describe('Browse Apps > Demo (Browse-DemoApps-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);

    await page.evaluate(() => (window as any).grok.shell.closeAll()).catch(() => undefined);
    await page.waitForTimeout(500);
    await ensureBrowsePanelOpen(page);
  });

  for (const demo of DEMOS) {
    const heavy = demo.tags?.includes('heavy');
    const id = demo.path.slice(2).join(' / ');
    const refNote = demo.ref ? ` [${demo.ref}]` : '';

    test(`Browse-DemoApps — ${id}${refNote}`, async ({ page }) => {
      test.setTimeout(heavy ? 150_000 : 60_000);

      if (demo.ref) test.info().annotations.push({ type: 'issue', description: demo.ref });
      const sink = watchErrors(page, demo.envIgnore ? [demo.envIgnore] : []);

      await openDemo(page, demo);

      await waitForDemoView(page, demo, heavy ? 90_000 : 30_000);

      await pokeView(page);

      await ensureContextPanelOpen(page, true);
      await expect(page.locator(CONTEXT_PANEL), 'Context Panel should be visible')
        .toBeVisible({ timeout: 10_000 });
      await page.waitForTimeout(1200);

      await expectNoErrors(page, sink);
    });
  }
});
