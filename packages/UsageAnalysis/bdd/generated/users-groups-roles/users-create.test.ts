/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-create.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, expand, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, dialogCloses, galleryCountLower, rememberGalleryCount, switchView, userStatusOnServer, usersOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
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
      await session.step(30, "And user types \"opavlenko+{time}c@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, session.text("opavlenko+{time}c@datagrok.ai"), el("Email input in \"Create new user\" dialog")));
      await session.step(31, "And user types \"opavlenko{time}c\" into Login input in \"Create new user\" dialog", () => typeInto(page, session.text("opavlenko{time}c"), el("Login input in \"Create new user\" dialog")));
      await session.step(32, "And user types \"Olesia\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "Olesia", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(33, "And user types \"BDD {time}c\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, session.text("BDD {time}c"), el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(34, "Then OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(35, "When user clicks on OK button in \"Create new user\" dialog", () => clickOn(page, el("OK button in \"Create new user\" dialog")));
      await session.step(36, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(37, "And the \"Olesia BDD {time}c\" view should be current", () => viewIsCurrent(page, session.text("Olesia BDD {time}c")));
      await session.step(38, "And 0 users with login \"opavlenko{time}c\" should be on the server", () => usersOnServer(page, 0, session.text("opavlenko{time}c")));
      await session.step(39, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(40, "Then 1 user with login \"opavlenko{time}c\" should be on the server", () => usersOnServer(page, 1, session.text("opavlenko{time}c")));
      await session.step(41, "And the user \"opavlenko{time}c\" should be active on the server", () => userStatusOnServer(page, session.text("opavlenko{time}c"), "active"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user is found in the Users view (Users-05)", async () => {
      await session.step(46, "When user switches to the \"Users\" view", () => switchView(page, "Users"));
      await session.step(47, "And user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(48, "And user types \"opavlenko{time}c\" into gallery search", () => typeInto(page, session.text("opavlenko{time}c"), el("gallery search")));
      await session.step(49, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(50, "And \"Olesia BDD {time}c\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"Olesia BDD {time}c\" link in gallery")), "visible"));
      await session.step(51, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A service user is saved by OK and shown its token (Users-07)", async () => {
      await session.step(56, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(57, "And user picks \"Service User...\" from the open menu", () => pickFromOpenMenu(page, "Service User..."));
      await session.step(58, "Then \"Create new service user\" dialog should be visible", () => shouldBe(page, el("\"Create new service user\" dialog"), "visible"));
      await session.step(59, "And OK button in \"Create new service user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "disabled"));
      await session.step(60, "When user types \"opavlenko-svc{time}c\" into Login input in \"Create new service user\" dialog", () => typeInto(page, session.text("opavlenko-svc{time}c"), el("Login input in \"Create new service user\" dialog")));
      await session.step(61, "Then OK button in \"Create new service user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new service user\" dialog"), "enabled"));
      await session.step(62, "When user clicks on OK button in \"Create new service user\" dialog", () => clickOn(page, el("OK button in \"Create new service user\" dialog")));
      await session.step(63, "Then the \"Create new service user\" dialog should close", () => dialogCloses(page, "Create new service user"));
      await session.step(64, "And \"API token\" dialog should be visible", () => shouldBe(page, el("\"API token\" dialog"), "visible"));
      await session.step(65, "And \"API token\" input in \"API token\" dialog should be visible", () => shouldBe(page, el("\"API token\" input in \"API token\" dialog"), "visible"));
      await session.step(66, "And 1 user with login \"opavlenko-svc{time}c\" should be on the server", () => usersOnServer(page, 1, session.text("opavlenko-svc{time}c")));
      await session.step(67, "When user clicks on CLOSE button in \"API token\" dialog", () => clickOn(page, el("CLOSE button in \"API token\" dialog")));
      await session.step(68, "Then the \"API token\" dialog should close", () => dialogCloses(page, "API token"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The service user is found in the Users view (Users-07)", async () => {
      await session.step(73, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(74, "And user types \"opavlenko-svc{time}c\" into gallery search", () => typeInto(page, session.text("opavlenko-svc{time}c"), el("gallery search")));
      await session.step(75, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(76, "And \"opavlenko-svc{time}c\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"opavlenko-svc{time}c\" link in gallery")), "visible"));
      await session.step(77, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(78, "Then no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
