/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/profile.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [core#1180]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dialogCloses, runningAccountSignedIn, signInAs, signInAsSelf, userNameOnServer, userOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A person edits their own name on their profile", () => {
  const session = feature(test, "features/general/profile.feature", import.meta.url);
  test("An edited name shows on the profile and is saved, and an empty or over-long first name is refused", {tag: ["@serial", "@realizes:core#1180"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
    await session.step(25, "When user signs in as \"bddprofile\"", () => signInAs(page, "bddprofile"));
    await session.step(26, "And user opens the address \"/u\"", () => openAddress(page, "/u"));
    await session.step(27, "Then \"Logout\" link should be visible", () => shouldBe(page, el("\"Logout\" link"), "visible"));
    await session.step(28, "When user clicks on second \"Edit\" icon", () => clickOn(page, el("second \"Edit\" icon")));
    await session.step(29, "Then \"Enter new name\" dialog should be visible", () => shouldBe(page, el("\"Enter new name\" dialog"), "visible"));
    await session.step(30, "And \"First name\" input in \"Enter new name\" dialog should have value \"bddprofile\"", () => shouldHaveValue(page, el("\"First name\" input in \"Enter new name\" dialog"), "bddprofile"));
    await session.step(31, "When user clears \"First name\" input in \"Enter new name\" dialog", () => clearField(page, el("\"First name\" input in \"Enter new name\" dialog")));
    await session.step(32, "Then OK button in \"Enter new name\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Enter new name\" dialog"), "disabled"));
    await session.step(33, "When user enters \"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL\" into \"First name\" input in \"Enter new name\" dialog", () => enterInto(page, "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL", el("\"First name\" input in \"Enter new name\" dialog")));
    await session.step(34, "Then \"First name\" input in \"Enter new name\" dialog should be invalid", () => shouldBe(page, el("\"First name\" input in \"Enter new name\" dialog"), "invalid"));
    await session.step(35, "And OK button in \"Enter new name\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Enter new name\" dialog"), "disabled"));
    await session.step(36, "When user enters \"Prof{time}\" into \"First name\" input in \"Enter new name\" dialog", () => enterInto(page, session.text("Prof{time}"), el("\"First name\" input in \"Enter new name\" dialog")));
    await session.step(37, "And user enters \"Ile{time}\" into \"Last name\" input in \"Enter new name\" dialog", () => enterInto(page, session.text("Ile{time}"), el("\"Last name\" input in \"Enter new name\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Enter new name\" dialog", () => clickOn(page, el("OK button in \"Enter new name\" dialog")));
    await session.step(39, "Then the \"Enter new name\" dialog should close", () => dialogCloses(page, "Enter new name"));
    await session.step(40, "And \"Prof{time} Ile{time}\" text should be visible", () => shouldBe(page, el(session.text("\"Prof{time} Ile{time}\" text")), "visible"));
    await session.step(41, "And the user \"bddprofile\" should have the name \"Prof{time} Ile{time}\" on the server", () => userNameOnServer(page, "bddprofile", session.text("Prof{time} Ile{time}")));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(43, "When user signs in as themselves again", () => signInAsSelf(page));
    await session.step(44, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
