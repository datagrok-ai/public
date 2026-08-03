import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions} from '../spec-login';

// Source scenario: General/profile-settings.md (automatable slice).
//
// The distinctive parts of the original scenario — uploading a profile photo
// in several image formats, rejecting a wrong-format file, and the
// "Change password..." dialog's client-side mismatch validation — are driven
// only through Dart-side UI paths (`user.setPicture`, the password Modal) with
// no public JS API surface, so they live in `profile-settings-ui.md`.
//
// What IS robustly automatable via the public API is the profile name edit and
// its persistence (the `ui.editableObject` name editor calls
// `dapi.users.save(user)`). This spec changes the current user's name through
// the same API the editor uses, verifies it persists across a re-fetch, and
// ALWAYS restores the original name in a finally block so the shared test
// account is left untouched.

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
    // --- Edit name and save (same path as the profile name editor).
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
    // friendlyName is derived from first + last names.
    expect(afterEdit.friendlyName).toContain(newFirst);
    expect(afterEdit.friendlyName).toContain(newLast);
  } finally {
    // --- Always restore the original name.
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
