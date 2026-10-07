/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/login-logout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {runningAccountSignedIn, sharingUserSignedIn, signInAsSelf, signInAsSharingUser} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Logging out, and the login form a signed-out page shows", () => {
  const session = feature(test, "features/general/login-logout.feature", import.meta.url);
  test("Logout leaves the shell for the login form, and an empty sign-in is refused there", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "When user signs in as the sharing user", () => signInAsSharingUser(page));
    await session.step(23, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
    await session.step(24, "When user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(25, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(26, "When user clicks on \"Logout\" link", () => clickOn(page, el("\"Logout\" link")));
    await session.step(27, "Then \"Login\" button should be visible", () => shouldBe(page, el("\"Login\" button"), "visible"));
    await session.step(28, "And \"Login failed\" text should be absent", () => shouldBe(page, el("\"Login failed\" text"), "absent"));
    await session.step(29, "And browse tab should be absent", () => shouldBe(page, el("browse tab"), "absent"));
    await session.step(30, "When user clicks on \"Login\" button", () => clickOn(page, el("\"Login\" button")));
    await session.step(31, "Then \"Login failed\" text should be visible", () => shouldBe(page, el("\"Login failed\" text"), "visible"));
    await session.step(32, "And \"Login\" button should be visible", () => shouldBe(page, el("\"Login\" button"), "visible"));
    await session.step(33, "And browse tab should be absent", () => shouldBe(page, el("browse tab"), "absent"));
    await session.step(34, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(35, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(36, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
  });
});
