/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, expand, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountHigher, galleryCountLower, rememberGalleryCount, switchView, userStatusOnServer, usersOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Creating users", () => {
  const session = feature(test, "features/users-groups-roles/users-create.feature", import.meta.url);
  test("Creating users", {tag: ["@journey", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(24, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await run.scenario("A new user is saved from the profile that OK opens (Users-05)", async () => {
      await session.step(27, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(28, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(29, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(30, "And user types \"bdd{time}@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd{time}@datagrok.ai"), el("Email input in \"Create new user\" dialog")));
      await session.step(31, "And user types \"bdd{time}\" into Login input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd{time}"), el("Login input in \"Create new user\" dialog")));
      await session.step(32, "And user types \"BDD\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "BDD", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(33, "And user types \"User {time}\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, session.text("User {time}"), el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(34, "Then OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(35, "When user clicks on OK button in \"Create new user\" dialog", () => clickOn(page, el("OK button in \"Create new user\" dialog")));
      await session.step(36, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(37, "And the \"BDD User {time}\" view should be current", () => viewIsCurrent(page, session.text("BDD User {time}")));
      await session.step(38, "And 0 users with login \"bdd{time}\" should be on the server", () => usersOnServer(page, 0, session.text("bdd{time}")));
      await session.step(39, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(40, "Then 1 user with login \"bdd{time}\" should be on the server", () => usersOnServer(page, 1, session.text("bdd{time}")));
      await session.step(41, "And the user \"bdd{time}\" should be active on the server", () => userStatusOnServer(page, session.text("bdd{time}"), "active"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user is found in the Users view (Users-05)", async () => {
      await session.step(46, "When user switches to the \"Users\" view", () => switchView(page, "Users"));
      await session.step(47, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(48, "And user types \"bdd{time}\" into gallery search", () => typeInto(page, session.text("bdd{time}"), el("gallery search")));
      await session.step(49, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(50, "And \"BDD User {time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD User {time}\" link in gallery")), "visible"));
      await session.step(51, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(52, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(53, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
      await session.step(55, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A service user is saved by OK and shown its token (Users-07)", async () => {
      await session.step(58, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(59, "And user picks \"Service User...\" from the open menu", () => pickFromOpenMenu(page, "Service User..."));
      await session.step(60, "Then \"Create new service user\" dialog should be visible", () => shouldBe(page, el("\"Create new service user\" dialog"), "visible"));
      await session.step(61, "And OK button in \"Create new service user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "disabled"));
      await session.step(62, "When user types \"bdd-svc{time}\" into Login input in \"Create new service user\" dialog", () => typeInto(page, session.text("bdd-svc{time}"), el("Login input in \"Create new service user\" dialog")));
      await session.step(63, "Then OK button in \"Create new service user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "enabled"));
      await session.step(64, "When user clicks on OK button in \"Create new service user\" dialog", () => clickOn(page, el("OK button in \"Create new service user\" dialog")));
      await session.step(65, "Then the \"Create new service user\" dialog should close", () => dialogCloses(page, "Create new service user"));
      await session.step(66, "And \"API token\" dialog should be visible", () => shouldBe(page, el("\"API token\" dialog"), "visible"));
      await session.step(67, "And \"API token\" input in \"API token\" dialog should be visible", () => shouldBe(page, el("\"API token\" input in \"API token\" dialog"), "visible"));
      await session.step(68, "And 1 user with login \"bdd-svc{time}\" should be on the server", () => usersOnServer(page, 1, session.text("bdd-svc{time}")));
      await session.step(69, "When user clicks on CLOSE button in \"API token\" dialog", () => clickOn(page, el("CLOSE button in \"API token\" dialog")));
      await session.step(70, "Then the \"API token\" dialog should close", () => dialogCloses(page, "API token"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
      await session.step(72, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The service user is found in the Users view (Users-07)", async () => {
      await session.step(75, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(76, "And user types \"bdd-svc{time}\" into gallery search", () => typeInto(page, session.text("bdd-svc{time}"), el("gallery search")));
      await session.step(77, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(78, "And \"bdd-svc{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bdd-svc{time}\" link in gallery")), "visible"));
      await session.step(79, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(80, "Then no errors should have been logged", () => noErrors(page));
      await session.step(81, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
