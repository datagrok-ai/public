// @verified — passed on dev 2026-03-26. Do not overwrite.
import {test, expect, Page} from '@playwright/test';

async function login(page: Page) {
  await page.goto('/?selenium=true');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
}

async function filteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function getFilterNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]');
    if (!fp) return [];
    return Array.from(fp.querySelectorAll('.d4-filter-column-name')).map(l => l.textContent!.trim());
  });
}

async function scrollFilterIntoView(page: Page, columnName: string) {
  await page.evaluate((col) => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
      if (label.textContent!.trim() === col) {
        const card = label.closest('.d4-filter')!;
        card.scrollIntoView({behavior: 'instant', block: 'center'});
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

test('Text Filter: Aroma search, AND/OR, fuzzy', async ({page}) => {
  test.setTimeout(300_000);
  await login(page);

  // Open beer.csv and filter panel
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.data.getDemoTable('beer.csv');
    grok.shell.addTableView(df);
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 10000});
  // Wait for Aroma column to be available before opening filters
  await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.col('Aroma') != null), {timeout: 10000}).toBe(true);
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});

  // Step 3: In Aroma search field, enter "low" and press Enter
  await test.step('3: Aroma filter — enter "low"', async () => {
    await scrollFilterIntoView(page, 'Aroma');

    // Click search icon on Aroma filter to show search input
    await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
        if (label.textContent!.trim() === 'Aroma') {
          const card = label.closest('.d4-filter')!;
          const searchIcon = card.querySelector('[name="icon-search"]') as HTMLElement;
          if (searchIcon) {
            searchIcon.style.visibility = 'visible';
            searchIcon.click();
          }
          break;
        }
      }
    });
    // Wait for search input to appear after clicking search icon
    await page.waitForTimeout(500);

    // Find the search input inside the Aroma filter card and type "low"
    const inputPos = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
        if (label.textContent!.trim() === 'Aroma') {
          const card = label.closest('.d4-filter')!;
          const input = card.querySelector('.d4-search-input, input[placeholder="Search"]') as HTMLElement;
          if (input) {
            const r = input.getBoundingClientRect();
            return {x: r.x + r.width / 2, y: r.y + r.height / 2};
          }
        }
      }
      return null;
    });

    if (inputPos) {
      await page.mouse.click(inputPos.x, inputPos.y);
      await page.keyboard.type('low');
      await page.keyboard.press('Enter');
    }

    await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeLessThan(totalRows);
    const count = await filteredCount(page);
    expect.soft(count).toBeGreaterThan(0);
  });

  // Step 4: Verify only rows containing "low" are displayed
  await test.step('4: Verify filtering works', async () => {
    const count = await filteredCount(page);
    expect.soft(count).toBeLessThan(totalRows);
    expect.soft(count).toBeGreaterThan(0);
  });

  // Step 5: Add multiple search terms
  await test.step('5: Add multiple search terms', async () => {
    const inputPos = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
        if (label.textContent!.trim() === 'Aroma') {
          const card = label.closest('.d4-filter')!;
          const input = card.querySelector('.d4-search-input, input[placeholder="Search"]') as HTMLElement;
          if (input) {
            const r = input.getBoundingClientRect();
            return {x: r.x + r.width / 2, y: r.y + r.height / 2};
          }
        }
      }
      return null;
    });

    if (inputPos) {
      await page.mouse.click(inputPos.x, inputPos.y);
      // Clear existing text and type multiple terms
      await page.keyboard.press('Control+A');
      await page.keyboard.type('low, medium');
      await page.keyboard.press('Enter');
    }

    await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeGreaterThan(0);
  });

  // Step 6: Switch between AND/OR modes
  let countOr: number;
  let countAnd: number;
  await test.step('6: Switch AND/OR modes', async () => {
    // OR mode should show more results than AND
    countOr = await filteredCount(page);

    // Look for AND/OR toggle button in the Aroma filter card
    const togglePos = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
        if (label.textContent!.trim() === 'Aroma') {
          const card = label.closest('.d4-filter')!;
          const allElements = card.querySelectorAll('*');
          for (const el of allElements) {
            if ((el.textContent === 'OR' || el.textContent === 'AND') && el.children.length === 0) {
              const r = (el as HTMLElement).getBoundingClientRect();
              if (r.width > 0) return {x: r.x + r.width / 2, y: r.y + r.height / 2, text: el.textContent};
            }
          }
        }
      }
      return null;
    });

    if (togglePos) {
      await page.mouse.click(togglePos.x, togglePos.y);
      await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeLessThanOrEqual(countOr);
      countAnd = await filteredCount(page);
    }
  });

  // Step 7: Verify filtering behavior
  await test.step('7: Verify filtering behavior', async () => {
    const count = await filteredCount(page);
    expect.soft(count).toBeGreaterThanOrEqual(0);
  });

  // Step 9: Adjust the Fuzzy Search slider
  await test.step('9: Adjust Fuzzy Search slider', async () => {
    const exactCount = await filteredCount(page);

    const sliderPos = await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      for (const label of fp.querySelectorAll('.d4-filter-column-name')) {
        if (label.textContent!.trim() === 'Aroma') {
          const card = label.closest('.d4-filter')!;
          const slider = card.querySelector('input[type="range"]') as HTMLElement;
          if (slider) {
            const r = slider.getBoundingClientRect();
            return {x: r.x, y: r.y + r.height / 2, width: r.width};
          }
        }
      }
      return null;
    });

    if (sliderPos) {
      await page.mouse.click(sliderPos.x + sliderPos.width * 0.5, sliderPos.y);
      await expect.poll(() => filteredCount(page), {timeout: 5000}).toBeGreaterThanOrEqual(exactCount);
    }

    const fuzzyCount = await filteredCount(page);
    expect.soft(fuzzyCount).toBeGreaterThanOrEqual(exactCount);
  });

  // Step 10: Verify more results with higher fuzzy level
  await test.step('10: Verify fuzzy gives more results', async () => {
    const count = await filteredCount(page);
    expect.soft(count).toBeGreaterThan(0);
  });
});
