/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/user-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
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
import {clearField, clickOn, enterInto, expand, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, galleryCountLower, rememberGalleryCount, userOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("User edge cases: the running account cannot disable itself, an invitation to a registered email", () => {
  const session = feature(test, "features/users-groups-roles/user-edge-cases.feature", import.meta.url);
  test("The running account's own context menu offers no Disable...", {tag: ["@users", "@edge-cases", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(21, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await session.step(24, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
    await session.step(25, "When user remembers the gallery counter", () => rememberGalleryCount(page));
    await session.step(26, "And user types \"{login}\" into gallery search", () => typeInto(page, session.text("{login}"), el("gallery search")));
    await session.step(27, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
    await session.step(28, "When user opens the context menu of \"{login}\" link in gallery", () => openContextMenu(page, el(session.text("\"{login}\" link in gallery"))));
    await session.step(29, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
    await session.step(30, "And the open menu should not list \"Disable...\"", () => menuDoesNotList(page, "Disable..."));
    await session.step(31, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(32, "And user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(33, "Then no errors should have been logged", () => noErrors(page));
    await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Invite a Friend with the email of a registered user says so", {tag: ["@users", "@edge-cases", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(21, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await session.step(37, "Given a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
    await session.step(38, "When user clicks on New button", () => clickOn(page, el("New button")));
    await session.step(39, "And user picks \"Invite a Friend...\" from the open menu", () => pickFromOpenMenu(page, "Invite a Friend..."));
    await session.step(40, "Then \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("\"Invite a Friend\" dialog"), "visible"));
    await session.step(41, "When user enters \"bddprofile@datagrok.ai\" into Email input in \"Invite a Friend\" dialog", () => enterInto(page, "bddprofile@datagrok.ai", el("Email input in \"Invite a Friend\" dialog")));
    await session.step(42, "Then \"Invite a Friend\" dialog should contain text \"User is already registered\"", () => shouldContainText(page, el("\"Invite a Friend\" dialog"), "User is already registered"));
    await session.step(43, "When user clicks on CANCEL button in \"Invite a Friend\" dialog", () => clickOn(page, el("CANCEL button in \"Invite a Friend\" dialog")));
    await session.step(44, "Then the \"Invite a Friend\" dialog should close", () => dialogCloses(page, "Invite a Friend"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
