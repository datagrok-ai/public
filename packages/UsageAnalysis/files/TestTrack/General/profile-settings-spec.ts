import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions} from '../spec-login';

test.use(specTestOptions);

test('Profile settings — name edit persists', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  const original = await page.evaluate(async () => {
    const grok = (window as any).grok;
    const u = await grok.dapi.users.current();
    return {id: u.id, firstName: u.firstName, lastName: u.lastName, friendlyName: u.friendlyName};
  });

  const newFirst = 'TtFirst' + Date.now().toString().slice(-6);
  const newLast = 'TtLast';

  try {

    const afterEdit = await page.evaluate(async ({first, last}) => {
      const grok = (window as any).grok;
      const u = await grok.dapi.users.current();
      u.firstName = first;
      u.lastName = last;
      await grok.dapi.users.save(u);
      const reloaded = await grok.dapi.users.find(u.id);
      return {firstName: reloaded.firstName, lastName: reloaded.lastName, friendlyName: reloaded.friendlyName};
    }, {first: newFirst, last: newLast});

    expect(afterEdit.firstName).toBe(newFirst);
    expect(afterEdit.lastName).toBe(newLast);

    expect(afterEdit.friendlyName).toContain(newFirst);
    expect(afterEdit.friendlyName).toContain(newLast);
  } finally {

    const restored = await page.evaluate(async (orig) => {
      const grok = (window as any).grok;
      const u = await grok.dapi.users.current();
      u.firstName = orig.firstName;
      u.lastName = orig.lastName;
      await grok.dapi.users.save(u);
      const reloaded = await grok.dapi.users.find(u.id);
      return {firstName: reloaded.firstName, lastName: reloaded.lastName};
    }, original);
    expect(restored.firstName).toBe(original.firstName);
    expect(restored.lastName).toBe(original.lastName);
  }
});
