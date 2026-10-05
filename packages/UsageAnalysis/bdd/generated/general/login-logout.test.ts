/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/login-logout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {runningAccountSignedIn, sharingUserSignedIn, signInAsSelf, signInAsSharingUser} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Logging out, and the login form a signed-out page shows", () => {
  const session = feature(test, "features/general/login-logout.feature", import.meta.url);
  test("Logout leaves the shell for the login form", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(27, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(28, "When user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(29, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(30, "When user clicks on \"Logout\" link", () => clickOn(page, el("\"Logout\" link")));
    await session.step(31, "Then \"Login\" button should be visible", () => shouldBe(page, el("\"Login\" button"), "visible"));
    await session.step(32, "And \"Login failed\" text should be absent", () => shouldBe(page, el("\"Login failed\" text"), "absent"));
    await session.step(33, "And browse tab should be absent", () => shouldBe(page, el("browse tab"), "absent"));
    await session.step(34, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(35, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(36, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
  });
  test("Login with both fields empty fails and keeps the form", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(39, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(40, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(41, "And user clicks on \"Logout\" link", () => clickOn(page, el("\"Logout\" link")));
    await session.step(42, "Then \"Login\" button should be visible", () => shouldBe(page, el("\"Login\" button"), "visible"));
    await session.step(43, "When user clicks on \"Login\" button", () => clickOn(page, el("\"Login\" button")));
    await session.step(44, "Then \"Login failed\" text should be visible", () => shouldBe(page, el("\"Login failed\" text"), "visible"));
    await session.step(45, "And \"Login\" button should be visible", () => shouldBe(page, el("\"Login\" button"), "visible"));
    await session.step(46, "And browse tab should be absent", () => shouldBe(page, el("browse tab"), "absent"));
    await session.step(47, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(48, "Then the running account should be signed in", () => runningAccountSignedIn(page));
  });
  test("A reload keeps the account signed in", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(52, "And user reloads the page", () => reloadPage(page));
    await session.step(53, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(54, "When user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(55, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(56, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
    await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(58, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(59, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
