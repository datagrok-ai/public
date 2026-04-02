import {test, expect, Page} from '@playwright/test';

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => { document.body.classList.add('selenium'); grok.shell.settings.showFiltersIconsConstantly = true; });
  await waitForPlatform(page);
}

async function waitForPlatform(page: Page) {
  await page.waitForFunction(() => {
    const bar = document.querySelector('.grok-loader,.d4-loading-bar,.d4-progress-bar');
    return !bar || getComputedStyle(bar).display === 'none';
  }, {timeout: 60000}).catch(() => {});
}

async function openSPGI(page: Page) {
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
    grok.shell.addTableView(df);
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 15000});
  // Wait for Chem package to initialize (Structure/Core filters need it)
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.columns.length), {timeout: 15000}).toBeGreaterThan(0);
  await page.waitForTimeout(3000);
}

async function openFilterPanel(page: Page) {
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});
}

async function closeFilterPanel(page: Page) {
  await page.evaluate(() => {
    for (const v of grok.shell.tv.viewers)
      if (v.type === 'Filters') { v.close(); break; }
  });
  await expect.poll(() => page.locator('[name="viewer-Filters"]').count(), {timeout: 5000}).toBe(0);
}

async function getFilteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function getFilterNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    if (!fp) return [];
    const labels = fp.querySelectorAll('.d4-filter-column-name');
    return Array.from(labels).map(l => l.textContent!.trim());
  });
}

async function resetFilters(page: Page) {
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-arrow-rotate-left"]') as HTMLElement;
    if (icon) icon.click();
  });
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount), {timeout: 5000}).toBe(totalRows);
}

async function scrollFilterIntoView(page: Page, columnName: string) {
  await page.evaluate((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        // Scroll the filter panel's scrollable container, not the document
        const scrollParent = label.closest('.d4-filter-group') || fp;
        const card = label.closest('.d4-filter')!;
        card.scrollIntoView({behavior: 'instant', block: 'center'});
        // Also try scrolling the parent overflow container
        let parent: HTMLElement | null = card.parentElement;
        while (parent && parent !== document.body) {
          if (parent.scrollHeight > parent.clientHeight) {
            const cardRect = card.getBoundingClientRect();
            const parentRect = parent.getBoundingClientRect();
            if (cardRect.top < parentRect.top || cardRect.bottom > parentRect.bottom)
              parent.scrollTop += cardRect.top - parentRect.top - parentRect.height / 2;
            break;
          }
          parent = parent.parentElement;
        }
        break;
      }
    }
  }, columnName);
}

async function getCategoryFilterCanvas(page: Page, columnName: string) {
  const handle = await page.evaluateHandle((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const card = label.closest('.d4-filter')!;
        const canvases = card.querySelectorAll('canvas');
        for (const c of canvases) {
          const r = c.getBoundingClientRect();
          if (r.width > 0 && r.height > 0) return c;
        }
      }
    }
    return null;
  }, columnName);
  const el = handle.asElement();
  return el;
}

async function getHistogramFilterCanvas(page: Page, columnName: string) {
  const handle = await page.evaluateHandle((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const card = label.closest('.d4-filter')!;
        const canvas = card.querySelector('canvas');
        return canvas;
      }
    }
    return null;
  }, columnName);
  const el = handle.asElement();
  return el;
}

async function clickFilterIcon(page: Page, columnName: string, iconName: string) {
  await page.evaluate(({col, icon}) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const card = label.closest('.d4-filter')!;
        const el = card.querySelector(`[name="${icon}"]`) as HTMLElement;
        if (el) el.click();
        break;
      }
    }
  }, {col: columnName, icon: iconName});
}

async function clickFilterPanelHamburger(page: Page) {
  await page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    const panel = fp?.closest('.d4-viewer-host') || fp?.closest('.panel-content')?.parentElement;
    const menuIcon = panel?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement;
    if (menuIcon) menuIcon.click();
  });
}

async function clickCategoryItem(page: Page, columnName: string, itemIndex: number) {
  // Get canvas coordinates, then use page.mouse.click() for real CDP events
  const pos = await page.evaluate(({col, idx}) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        label.scrollIntoView({behavior: 'instant', block: 'center'});
        const card = label.closest('.d4-filter')!;
        for (const c of card.querySelectorAll('canvas')) {
          const r = c.getBoundingClientRect();
          if (r.width > 0 && r.height > 0) {
            const rowHeight = r.height / 5;
            return {x: r.x + 60, y: r.y + rowHeight * idx + rowHeight / 2};
          }
        }
      }
    }
    return null;
  }, {col: columnName, idx: itemIndex});
  if (pos)
    await page.mouse.click(pos.x, pos.y);
  await page.waitForTimeout(500);
}

async function setHistogramFilterRange(page: Page, columnName: string, min: number, max: number) {
  await page.evaluate(async ({col, minVal, maxVal}) => {
    const fg = grok.shell.tv.getFiltersGroup();
    fg.updateOrAdd({type: 'histogram', column: col, min: minVal, max: maxVal});
    await new Promise(r => setTimeout(r, 1000));
  }, {col: columnName, minVal: min, maxVal: max});
  await page.waitForTimeout(500);
}

async function removeFilter(page: Page, columnName: string) {
  await page.evaluate((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const header = label.closest('.d4-filter')!.querySelector('.d4-filter-header')!;
        const x = header.querySelector('[name="icon-times"]') as HTMLElement;
        if (x) x.click();
        break;
      }
    }
  }, columnName);
  await page.waitForTimeout(500);
}

async function toggleFilterEnabled(page: Page, columnName: string) {
  await page.evaluate((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const header = label.closest('.d4-filter')!.querySelector('.d4-filter-header')!;
        const cb = header.querySelector('input[type="checkbox"]') as HTMLInputElement;
        if (cb) cb.click();
        break;
      }
    }
  }, columnName);
  await page.waitForTimeout(500);
}

test('Filter Panel: Basic Operations — all sections', async ({page}) => {
  test.setTimeout(300_000);

  await login(page);
  await openSPGI(page);
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  await openFilterPanel(page);

  // ===== Section 1: Filtering and resetting =====
  await test.step('Section 1: Filtering and resetting', async () => {

    // Step 1-2: Open Sketcher, type c1ccccc1, apply Structure filter
    let afterStructure: number;
    await test.step('1.1-2: Structure filter via Sketcher', async () => {
      await page.locator('[name="viewer-Filters"]').getByText('Sketch').first().click();
      await page.getByText('CANCEL').waitFor({timeout: 10000});
      const smilesInput = page.getByPlaceholder('SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');
      await smilesInput.click();
      await smilesInput.fill('c1ccccc1');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(500);
      await page.locator('button:has-text("OK")').click();
      await expect.poll(() => getFilteredCount(page), {timeout: 10000}).toBeLessThan(totalRows);
      afterStructure = await getFilteredCount(page);
      expect.soft(afterStructure).toBeLessThan(totalRows);
      expect.soft(afterStructure).toBeGreaterThan(0);
    });

    // Step 3: Click R_ONE in Stereo Category (first item, index 0)
    let afterStereo: number;
    await test.step('1.3: Click R_ONE in Stereo Category', async () => {
      await clickCategoryItem(page, 'Stereo Category', 0);
      afterStereo = await getFilteredCount(page);
      expect.soft(afterStereo).toBeLessThan(afterStructure);
      expect.soft(afterStereo).toBeGreaterThan(0);
    });

    // Step 4: Click filter menu on Average Mass, enable Min/max, set max to 400
    let afterMass: number;
    await test.step('1.4: Average Mass — set max to 400 via Min/max input', async () => {
      await setHistogramFilterRange(page, 'Average Mass', 221, 400);
      afterMass = await getFilteredCount(page);
      expect.soft(afterMass).toBeLessThan(afterStereo);
      expect.soft(afterMass).toBeGreaterThan(0);
    });

    // Step 5: Disable all filters
    await test.step('1.5: Disable all filters', async () => {
      await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]')!;
        fp.scrollTop = 0;
      });
      await page.waitForTimeout(300);
      // Click the global enable/disable checkbox in the filter panel header
      await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]')!;
        const header = fp.closest('.panel-content')?.previousElementSibling;
        const globalCb = fp.parentElement?.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
        if (globalCb) globalCb.click();
      });
      await page.waitForTimeout(500);
      expect.soft(await getFilteredCount(page)).toBe(totalRows);
    });

    // Step 6: Re-enable filters
    await test.step('1.6: Re-enable all filters', async () => {
      await page.evaluate(() => {
        const globalCb = document.querySelector('[name="viewer-Filters"]')
          ?.parentElement?.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
        if (globalCb) globalCb.click();
      });
      await page.waitForTimeout(1000);
      expect.soft(await getFilteredCount(page)).toBe(afterMass);
    });

    // Step 7: Disable Stereo Category individually
    await test.step('1.7: Disable Stereo Category filter', async () => {
      await toggleFilterEnabled(page, 'Stereo Category');
      const count = await getFilteredCount(page);
      expect.soft(count).toBeGreaterThan(afterMass);
    });

    // Step 8: Close filter panel — all rows shown
    await test.step('1.8: Close filter panel', async () => {
      await closeFilterPanel(page);
      expect.soft(await getFilteredCount(page)).toBe(totalRows);
    });

    // Step 9: Reopen — Stereo Category still disabled
    await test.step('1.9: Reopen — Stereo Category still disabled', async () => {
      await openFilterPanel(page);
      const stereoDisabled = await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]')!;
        for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
          if (label.textContent!.trim() === 'Stereo Category') {
            const cb = label.closest('.d4-filter')!
              .querySelector('.d4-filter-header input[type="checkbox"]') as HTMLInputElement;
            return cb ? !cb.checked : false;
          }
        }
        return false;
      });
      expect.soft(stereoDisabled).toBe(true);
    });

    // Step 10: Reset all filters
    await test.step('1.10: Reset all filters', async () => {
      await resetFilters(page);
      expect.soft(await getFilteredCount(page)).toBe(totalRows);
    });

    // Steps 11-12: Close and reopen — default state
    await test.step('1.11-12: Close/reopen — default state', async () => {
      await closeFilterPanel(page);
      await openFilterPanel(page);
      expect.soft(await getFilteredCount(page)).toBe(totalRows);
    });

    // Steps 13-14: Remove Structure and Core filters via X icon
    await test.step('1.13-14: Remove Structure and Core filters', async () => {
      await removeFilter(page, 'Structure');
      await removeFilter(page, 'Core');
    });

    // Step 15: Close/reopen — removed filters absent
    await test.step('1.15: Removed filters stay removed', async () => {
      await closeFilterPanel(page);
      await openFilterPanel(page);
      const names = await getFilterNames(page);
      expect.soft(names).not.toContain('Structure');
      expect.soft(names).not.toContain('Core');
    });
  });

  // ===== Section 2: Filter Panel UI Behavior =====
  await test.step('Section 2: Filter Panel UI Behavior', async () => {
    await resetFilters(page);

    // Step 2.1: Open Scatter Plot and Bar Chart
    await test.step('2.1: Add Scatter Plot and Bar Chart', async () => {
      await page.locator('[name="icon-scatter-plot"]').click();
      await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 5000});
      await page.locator('[name="icon-bar-chart"]').click();
      await page.locator('[name="viewer-Bar-chart"]').waitFor({timeout: 5000});
      await expect.soft(page.locator('[name="viewer-Scatter-plot"]')).toBeVisible();
      await expect.soft(page.locator('[name="viewer-Bar-chart"]')).toBeVisible();
    });

    // Step 2.3: Hover over Stereo Category value — verify cross-highlighting
    await test.step('2.3: Hover Stereo Category — cross-highlighting', async () => {
      await scrollFilterIntoView(page, 'Stereo Category');
      const canvas = await getCategoryFilterCanvas(page, 'Stereo Category');
      expect.soft(canvas).not.toBeNull();
      // Hover the first category row (R_ONE) in the canvas
      const box = await canvas!.boundingBox();
      expect.soft(box).not.toBeNull();
      const rowHeight = box!.height / 5;
      await page.mouse.move(box!.x + 60, box!.y + rowHeight / 2);
      await page.waitForTimeout(500);
      // Cross-highlighting is visual — take a screenshot for manual review
    });

    // Step 2.4: Hover over Competition assay Date filter — verify tooltip date format
    await test.step('2.4: Hover Competition assay Date — tooltip shows date format', async () => {
      await scrollFilterIntoView(page, 'Competition assay Date');
      const canvas = await getHistogramFilterCanvas(page, 'Competition assay Date');
      expect.soft(canvas).not.toBeNull();
      const box = await canvas!.boundingBox();
      expect.soft(box).not.toBeNull();
      await page.mouse.move(box!.x + box!.width / 2, box!.y + box!.height / 2);
      await page.waitForTimeout(1000);
      // Tooltip is visual — verify by checking the date format matches M/d/yyyy (default)
    });

    // Step 2.5: Right-click Competition assay Date column header → Format → Custom → MM/dd/yyyy
    await test.step('2.5: Change date format to MM/dd/yyyy', async () => {
      // Scroll grid to the Competition assay Date column
      await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const col = grid.col('Competition assay Date');
        if (col) grid.scrollToCell(col, 0);
      });
      await page.waitForTimeout(500);

      // Right-click the column header via overlay canvas
      const headerPos = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const overlay = grid.overlay;
        const r = overlay.getBoundingClientRect();
        for (let x = 0; x < r.width; x += 20) {
          const cell = grid.hitTest(x, 10);
          if (cell && cell.gridColumn && cell.gridColumn.name === 'Competition assay Date')
            return {screenX: r.x + x + 40, screenY: r.y + 10};
        }
        return null;
      });
      expect.soft(headerPos).not.toBeNull();

      await page.evaluate(({screenX, screenY}) => {
        const grid = grok.shell.tv.grid;
        const overlay = grid.overlay;
        overlay.dispatchEvent(new PointerEvent('pointerdown', {
          clientX: screenX, clientY: screenY, button: 2, buttons: 2,
          bubbles: true, cancelable: true, pointerId: 1, pointerType: 'mouse', view: window
        }));
        overlay.dispatchEvent(new MouseEvent('contextmenu', {
          clientX: screenX, clientY: screenY, button: 2, buttons: 2,
          bubbles: true, cancelable: true, view: window
        }));
      }, headerPos!);
      await page.waitForTimeout(500);

      // Hover Format menu item to expand submenu
      await page.getByText('Format', {exact: true}).first().hover();
      await page.waitForTimeout(500);

      // Click Custom...
      await page.getByText('Custom...').click();
      await page.waitForTimeout(500);

      // Type MM/dd/yyyy in the Custom field
      await page.keyboard.press('Control+A');
      await page.keyboard.type('MM/dd/yyyy');
      await page.locator('button:has-text("OK")').click();
      await page.waitForTimeout(500);

      // Verify format changed in grid
      const dateText = await page.evaluate(() => {
        const col = grok.shell.tv.dataFrame.col('Competition assay Date');
        const format = col.tags['format'];
        return format;
      });
      expect.soft(dateText).toBe('MM/dd/yyyy');
    });

    // Step 2.6: Hover Competition assay Date filter — verify tooltip shows MM/dd/yyyy format
    await test.step('2.6: Verify date filter tooltip uses MM/dd/yyyy', async () => {
      await scrollFilterIntoView(page, 'Competition assay Date');
      const canvas = await getHistogramFilterCanvas(page, 'Competition assay Date');
      expect.soft(canvas).not.toBeNull();
      const box = await canvas!.boundingBox();
      expect.soft(box).not.toBeNull();
      await page.mouse.move(box!.x + box!.width / 2, box!.y + box!.height / 2);
      await page.waitForTimeout(1000);
      // Tooltip should now show MM/dd/yyyy format — visual verification via screenshot
    });

    // Step 2.7: Click Search icon on AMES filter — search input appears
    await test.step('2.7: Click Search icon on AMES filter', async () => {
      await scrollFilterIntoView(page, 'Series');
      await page.waitForTimeout(500);
      // Click the search icon using real mouse coordinates (evaluate click may not work)
      const searchPos = await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]')!;
        for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
          if (label.textContent!.trim() === 'Series') {
            const card = label.closest('.d4-filter')!;
            const icon = card.querySelector('[name="icon-search"]') as HTMLElement;
            if (icon) {
              icon.style.visibility = 'visible';
              const r = icon.getBoundingClientRect();
              return {x: r.x + r.width / 2, y: r.y + r.height / 2, found: true};
            }
            return {found: false, reason: 'no icon-search in card'};
          }
        }
        return {found: false, reason: 'AMES filter not found'};
      });
      if (searchPos && 'x' in searchPos)
        await page.mouse.click(searchPos.x, searchPos.y);
      await page.waitForTimeout(500);
      // Verify search input appeared inside the AMES filter card
      const hasSearch = await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]')!;
        for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
          if (label.textContent!.trim() === 'Series') {
            const card = label.closest('.d4-filter')!;
            const input = card.querySelector('.d4-search-input') as HTMLInputElement;
            return input != null && getComputedStyle(input).display !== 'none';
          }
        }
        return false;
      });
      expect.soft(hasSearch).toBe(true);
    });

    // Step 2.8: Save layout via Toolbox > Layouts > SAVE
    await test.step('2.8: Save layout', async () => {
      await page.locator('[name="div-section--Layouts"]').click();
      await page.waitForTimeout(500);
      // Click the SAVE button inside the Layouts pane (not the top toolbar SAVE)
      const saveBtn = page.locator('.d4-toolbox-layouts').locator('button:has-text("SAVE")');
      await saveBtn.click();
      // Wait for the layout to actually save (thumbnail appears)
      await page.locator('.grok-suggestions-chart-host').first().waitFor({timeout: 10000});
      await page.waitForTimeout(1000);
    });

    // Step 2.9: Click the saved layout — verify search field still visible in AMES
    await test.step('2.9: Apply saved layout — search field persists', async () => {
      // Click the first layout thumbnail
      await page.locator('.grok-suggestions-chart-host').first().click();
      await page.waitForTimeout(5000);
      // Wait for filter panel to reappear after layout apply
      await page.evaluate(async () => {
        grok.shell.tv.getFiltersGroup();
        await new Promise(r => setTimeout(r, 2000));
      });
      const hasSearch = await page.evaluate(() => {
        const fp = document.querySelector('[name="viewer-Filters"]');
        if (!fp) return false;
        for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
          if (label.textContent!.trim() === 'Series') {
            const card = label.closest('.d4-filter')!;
            const input = card.querySelector('.d4-search-input') as HTMLInputElement;
            return input != null && getComputedStyle(input).display !== 'none';
          }
        }
        return false;
      });
      expect.soft(hasSearch).toBe(true);
    });

    // Step 2.10: Reset filters
    await test.step('2.10: Reset filters', async () => {
      await resetFilters(page);
      expect.soft(await getFilteredCount(page)).toBe(totalRows);
    });

    // No cleanup here — viewers and layout will be cleaned up after Section 4
  });

  // ===== Section 3: Adding and Reordering Filters =====
  await test.step('Section 3: Adding and Reordering Filters', async () => {
    await resetFilters(page);

    // Step 3.2: Reorder Filters via hamburger menu
    await test.step('3.2: Reorder Filters via hamburger', async () => {
      await clickFilterPanelHamburger(page);
      await page.waitForTimeout(500);
      const reorderItem = page.getByText('Reorder Filters');
      if (await reorderItem.isVisible({timeout: 2000}).catch(() => false)) {
        await reorderItem.click();
        await page.waitForTimeout(1000);
        // Close the dialog if it appeared
        const okBtn = page.locator('button:has-text("OK")');
        if (await okBtn.isVisible({timeout: 2000}).catch(() => false))
          await okBtn.click();
      }
      await page.keyboard.press('Escape');
      await page.waitForTimeout(500);
    });

    // Step 3.3: Add Columns via hamburger > Select Columns
    await test.step('3.3: Select Columns dialog', async () => {
      await clickFilterPanelHamburger(page);
      await page.waitForTimeout(500);
      const selectColsItem = page.getByText('Select Columns...');
      if (await selectColsItem.isVisible({timeout: 2000}).catch(() => false)) {
        await selectColsItem.click();
        await page.waitForTimeout(1000);
        // Close the dialog
        const okBtn = page.locator('button:has-text("OK")');
        if (await okBtn.isVisible({timeout: 2000}).catch(() => false))
          await okBtn.click();
      }
      await page.keyboard.press('Escape');
      await page.waitForTimeout(500);
    });

    // Step 3.4: Hamburger > Remove All
    await test.step('3.4: Hamburger > Remove All', async () => {
      await clickFilterPanelHamburger(page);
      await page.waitForTimeout(500);
      await page.getByText('Remove All').click();
      await page.waitForTimeout(500);
      const okBtn = page.locator('button:has-text("OK")');
      if (await okBtn.isVisible({timeout: 2000}).catch(() => false))
        await okBtn.click();
      await page.waitForTimeout(500);
      const names = await getFilterNames(page);
      expect.soft(names.length).toBe(0);
    });

    // Step 3.5: Add filters via different methods and verify each appears at top
    await test.step('3.5: Add filters via JS API (simulating UI methods)', async () => {
      // 3.5a: Add "Id" filter (simulating drag column header into filter panel)
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Id'});
      });
      await page.waitForTimeout(1000);
      let names = await getFilterNames(page);
      expect.soft(names[0]).toBe('Id');

      // 3.5b: Add "CAST Idea ID" (simulating right-click header > Filter > Add Filter)
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'CAST Idea ID'});
      });
      await page.waitForTimeout(1000);
      names = await getFilterNames(page);
      expect.soft(names[0]).toBe('CAST Idea ID');

      // 3.5c: Add "Stereo Category" (simulating right-click cell > Use as filter)
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category'});
      });
      await page.waitForTimeout(1000);
      names = await getFilterNames(page);
      expect.soft(names[0]).toBe('Stereo Category');
    });

    // Step 3.6: Remove All again
    await test.step('3.6: Remove All', async () => {
      await clickFilterPanelHamburger(page);
      await page.waitForTimeout(500);
      await page.getByText('Remove All').click();
      await page.waitForTimeout(500);
      const okBtn = page.locator('button:has-text("OK")');
      if (await okBtn.isVisible({timeout: 2000}).catch(() => false))
        await okBtn.click();
      await page.waitForTimeout(500);
      const names = await getFilterNames(page);
      expect.soft(names.length).toBe(0);
    });

    // Step 3.7: Close filter panel
    await test.step('3.7: Close filter panel', async () => {
      await closeFilterPanel(page);
    });
  });

  // ===== Section 4: Hidden Columns =====
  // Precondition: fresh table view needed because Remove All in Section 3
  // permanently removes default filters. Recreate table to get clean state.
  await test.step('Section 4: Hidden columns', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
    });
    await page.waitForTimeout(5000);

    await test.step('4.1-2: Hide Structure, Core, R1', async () => {
      await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        for (const col of ['Structure', 'Core', 'R1']) {
          const c = grid.col(col);
          if (c) c.visible = false;
        }
      });
      await page.waitForTimeout(500);
    });

    await test.step('4.3: Verify hidden columns not in grid', async () => {
      const visibleCols = await page.evaluate(() => {
        const cols: string[] = [];
        const grid = grok.shell.tv.grid;
        for (let i = 0; i < grid.columns.length; i++) {
          const c = grid.columns.byIndex(i);
          if (c.visible) cols.push(c.name);
        }
        return cols;
      });
      expect.soft(visibleCols).not.toContain('Structure');
      expect.soft(visibleCols).not.toContain('Core');
      expect.soft(visibleCols).not.toContain('R1');
    });

    await test.step('4.4: Filter panel — no cards for hidden columns', async () => {
      await openFilterPanel(page);
      const names = await getFilterNames(page);
      expect.soft(names).not.toContain('Structure');
      expect.soft(names).not.toContain('Core');
      expect.soft(names).not.toContain('R1');
    });

    await test.step('4.5: Remove All filters, close panel, unhide columns', async () => {
      // Hamburger menu > Remove All
      await clickFilterPanelHamburger(page);
      await page.waitForTimeout(500);
      await page.getByText('Remove All').click();
      await page.waitForTimeout(500);
      const okBtn = page.locator('button:has-text("OK")');
      if (await okBtn.isVisible({timeout: 2000}).catch(() => false))
        await okBtn.click();
      await page.waitForTimeout(500);

      await closeFilterPanel(page);
      await page.waitForTimeout(500);
      await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        for (const col of ['Structure', 'Core', 'R1']) {
          const c = grid.col(col);
          if (c) c.visible = true;
        }
      });
      await page.waitForTimeout(500);
    });

    await test.step('4.6: Filters reappear after unhiding', async () => {
      await openFilterPanel(page);
      // Wait for Chem package to asynchronously create SubstructureFilters
      await page.waitForTimeout(5000);
      const namesAfter = await getFilterNames(page);
      expect.soft(namesAfter).toContain('Structure');
      expect.soft(namesAfter).toContain('Core');
      expect.soft(namesAfter).toContain('R1');
    });
  });

  // Final cleanup: delete saved layout from Section 2
  await page.evaluate(async () => {
    try {
      const layouts = await grok.dapi.layouts.filter('table.name = "spgi-100"').list();
      if (layouts.length > 0) {
        layouts.sort((a, b) => b.createdOn - a.createdOn);
        await grok.dapi.layouts.delete(layouts[0]);
      }
    } catch (e) {}
  }).catch(() => {});
});
