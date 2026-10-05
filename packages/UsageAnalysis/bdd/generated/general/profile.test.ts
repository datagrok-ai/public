/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/profile.feature
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
import {clearField, clickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dialogCloses, runningAccountSignedIn, signInAs, signInAsSelf, userOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A person edits their own name on their profile", () => {
  const session = feature(test, "features/general/profile.feature", import.meta.url);
  test("An edited name is saved, survives a reload, and is put back", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
    await session.step(21, "When user signs in as \"bddprofile\"", () => signInAs(page, "bddprofile"));
    await session.step(22, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(23, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(24, "When user clicks on second \"Edit\" icon", () => clickOn(page, el("second \"Edit\" icon")));
    await session.step(25, "Then \"Enter new name\" dialog should be visible", () => shouldBe(page, el("\"Enter new name\" dialog"), "visible"));
    await session.step(26, "When user enters \"Prof{time}\" into \"First name\" input in \"Enter new name\" dialog", () => enterInto(page, session.text("Prof{time}"), el("\"First name\" input in \"Enter new name\" dialog")));
    await session.step(27, "And user enters \"Ile{time}\" into \"Last name\" input in \"Enter new name\" dialog", () => enterInto(page, session.text("Ile{time}"), el("\"Last name\" input in \"Enter new name\" dialog")));
    await session.step(28, "And user clicks on OK button in \"Enter new name\" dialog", () => clickOn(page, el("OK button in \"Enter new name\" dialog")));
    await session.step(29, "Then the \"Enter new name\" dialog should close", () => dialogCloses(page, "Enter new name"));
    await session.step(30, "When user reloads the page", () => reloadPage(page));
    await session.step(31, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(32, "Then \"Prof{time} Ile{time}\" text should be visible", () => shouldBe(page, el(session.text("\"Prof{time} Ile{time}\" text")), "visible"));
    await session.step(33, "When user clicks on second \"Edit\" icon", () => clickOn(page, el("second \"Edit\" icon")));
    await session.step(34, "Then \"Enter new name\" dialog should be visible", () => shouldBe(page, el("\"Enter new name\" dialog"), "visible"));
    await session.step(35, "When user enters \"bddprofile\" into \"First name\" input in \"Enter new name\" dialog", () => enterInto(page, "bddprofile", el("\"First name\" input in \"Enter new name\" dialog")));
    await session.step(36, "And user clears \"Last name\" input in \"Enter new name\" dialog", () => clearField(page, el("\"Last name\" input in \"Enter new name\" dialog")));
    await session.step(37, "And user clicks on OK button in \"Enter new name\" dialog", () => clickOn(page, el("OK button in \"Enter new name\" dialog")));
    await session.step(38, "Then the \"Enter new name\" dialog should close", () => dialogCloses(page, "Enter new name"));
    await session.step(39, "When user reloads the page", () => reloadPage(page));
    await session.step(40, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(41, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(42, "And \"Prof{time} Ile{time}\" text should be absent", () => shouldBe(page, el(session.text("\"Prof{time} Ile{time}\" text")), "absent"));
    await session.step(43, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(44, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
