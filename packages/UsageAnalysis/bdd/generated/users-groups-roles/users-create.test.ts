/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, close, expand, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountHigher, galleryCountLower, rememberGalleryCount, userOnServer, usersOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating users", () => {
  const session = feature(test, "features/users-groups-roles/users-create.feature", import.meta.url);
  test("Creating users", {tag: ["@journey", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And a user \"bddcreated\" is on the server", () => userOnServer(page, "bddcreated"));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(24, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await run.scenario("OK in the new user dialog opens the new user's profile and saves nobody (Users-05)", async () => {
      await session.step(27, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(28, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(29, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(30, "And user types \"bddunsaved@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, "bddunsaved@datagrok.ai", el("Email input in \"Create new user\" dialog")));
      await session.step(31, "And user types \"bddunsaved\" into Login input in \"Create new user\" dialog", () => typeInto(page, "bddunsaved", el("Login input in \"Create new user\" dialog")));
      await session.step(32, "And user types \"BDD\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "BDD", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(33, "And user types \"Unsaved\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, "Unsaved", el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(34, "Then OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(35, "When user clicks on OK button in \"Create new user\" dialog", () => clickOn(page, el("OK button in \"Create new user\" dialog")));
      await session.step(36, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(37, "And the \"BDD Unsaved\" view should be current", () => viewIsCurrent(page, "BDD Unsaved"));
      await session.step(38, "And 0 users with login \"bddunsaved\" should be on the server", () => usersOnServer(page, 0, "bddunsaved"));
      await session.step(39, "When user closes \"BDD Unsaved\" view", () => close(page, el("\"BDD Unsaved\" view")));
      await session.step(40, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(41, "And 0 users with login \"bddunsaved\" should be on the server", () => usersOnServer(page, 0, "bddunsaved"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user is found in the Users view (Users-05)", async () => {
      await session.step(46, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(47, "And user types \"bddcreated\" into gallery search", () => typeInto(page, "bddcreated", el("gallery search")));
      await session.step(48, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(49, "And \"bddcreated\" link in gallery should be visible", () => shouldBe(page, el("\"bddcreated\" link in gallery"), "visible"));
      await session.step(50, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(51, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(52, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new service user dialog asks for a login before OK (Users-07)", async () => {
      await session.step(57, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(58, "And user picks \"Service User...\" from the open menu", () => pickFromOpenMenu(page, "Service User..."));
      await session.step(59, "Then \"Create new service user\" dialog should be visible", () => shouldBe(page, el("\"Create new service user\" dialog"), "visible"));
      await session.step(60, "And OK button in \"Create new service user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "disabled"));
      await session.step(61, "When user types \"bdd-svc-unsaved\" into Login input in \"Create new service user\" dialog", () => typeInto(page, "bdd-svc-unsaved", el("Login input in \"Create new service user\" dialog")));
      await session.step(62, "Then OK button in \"Create new service user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "enabled"));
      await session.step(63, "When user clicks on CANCEL button in \"Create new service user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new service user\" dialog")));
      await session.step(64, "Then the \"Create new service user\" dialog should close", () => dialogCloses(page, "Create new service user"));
      await session.step(65, "And 0 users with login \"bdd-svc-unsaved\" should be on the server", () => usersOnServer(page, 0, "bdd-svc-unsaved"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
