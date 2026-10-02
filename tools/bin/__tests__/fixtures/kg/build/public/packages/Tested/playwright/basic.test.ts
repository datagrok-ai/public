import {test, expect} from '@playwright/test';

test.describe('Basic', () => {
  test.beforeEach(async ({page}) => {
    await page.goto('/');
  });

  test('~domains/bio sequence view opens', async ({page}) => {
    await expect(page.locator('.grok-view')).toBeVisible();
  });

  test.describe('Inner', () => {
    test('nested title', async ({page}) => {
      const s = "it's } tricky";
    });
  });
});

test.skip('later', async ({page}) => {});

test.describe.skip('Parked', () => {
  test('parked test', async ({page}) => {});
});

// test('commented', async ({page}) => {});
