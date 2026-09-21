/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/roles-view.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.roles]
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
import {clearField, clickOn, expand, followingShouldBe, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dialogCloses, galleryCountHigher, galleryCountLower, galleryMode, groupsOnServer, noGroupOnServer, rememberGalleryCount, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Roles view", () => {
  const session = feature(test, "features/users-groups-roles/roles-view.feature", import.meta.url);
  test("The Roles view", {tag: ["@journey", "@serial", "@roles", "@realizes:views.roles"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And no role named \"BDD-RV-Role-{time}\" is on the server", () => noGroupOnServer(page, session.text("BDD-RV-Role-{time}")));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(23, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(24, "And user clicks on \"Platform > Roles\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Roles\" tree node inside browse tree")));
    await run.scenario("The view opens from the Browse tree (Roles-01)", async () => {
      await session.step(27, "Then the \"Roles\" view should be current", () => viewIsCurrent(page, "Roles"));
      await session.step(28, "And the page address should contain \"/roles\"", () => urlShouldContain(page, "/roles"));
      await session.step(29, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(30, "And gallery counter should be visible", () => shouldBe(page, el("gallery counter"), "visible"));
      await session.step(31, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A role to look at", async () => {
      await session.step(36, "When user clicks on \"New Role...\" button", () => clickOn(page, el("\"New Role...\" button")));
      await session.step(37, "And user types \"BDD-RV-Role-{time}\" into Name input in \"Create New Role\" dialog", () => typeInto(page, session.text("BDD-RV-Role-{time}"), el("Name input in \"Create New Role\" dialog")));
      await session.step(38, "And user clicks on OK button in \"Create New Role\" dialog", () => clickOn(page, el("OK button in \"Create New Role\" dialog")));
      await session.step(39, "Then the \"Create New Role\" dialog should close", () => dialogCloses(page, "Create New Role"));
      await session.step(40, "And 1 role named \"BDD-RV-Role-{time}\" should be on the server", () => groupsOnServer(page, 1, session.text("BDD-RV-Role-{time}")));
      await session.step(41, "When user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(42, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(43, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(44, "Then no errors should have been logged", () => noErrors(page));
      await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The toolbar carries its controls (Roles-02)", async () => {
      await session.step(48, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"New Role...\" button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]));
      await session.step(55, "And \"New Role...\" button should be enabled", () => shouldBe(page, el("\"New Role...\" button"), "enabled"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
      await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The view-mode icons switch the gallery's render mode (Roles-07)", async () => {
      await session.step(60, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(61, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(62, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(63, "And \"Switch to card view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(64, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(65, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(66, "And \"Switch to grid view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to grid view\" icon inside gallery toolbar"), "selected"));
      await session.step(67, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(68, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(69, "And \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
      await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Searching by name finds the one role, and clearing brings the rest back (Roles-06)", async () => {
      await session.step(74, "When user types \"BDD-RV-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RV-Role-{time}"), el("gallery search")));
      await session.step(75, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(76, "And \"BDD-RV-Role-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-RV-Role-{time}\" link in gallery")), "visible"));
      await session.step(77, "And the page address should contain \"?q=BDD-RV-Role-{time}\"", () => urlShouldContain(page, session.text("?q=BDD-RV-Role-{time}")));
      await session.step(78, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(79, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(80, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(81, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(82, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(83, "Then no errors should have been logged", () => noErrors(page));
      await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A role's context menu (Roles-08)", async () => {
      await session.step(87, "When user types \"BDD-RV-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RV-Role-{time}"), el("gallery search")));
      await session.step(88, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(89, "When user opens the context menu of \"BDD-RV-Role-{time}\" link in gallery", () => openContextMenu(page, el(session.text("\"BDD-RV-Role-{time}\" link in gallery"))));
      await session.step(90, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(91, "And the open menu should list \"Delete\"", () => menuLists(page, "Delete"));
      await session.step(92, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(93, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(94, "Then no errors should have been logged", () => noErrors(page));
      await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting a role fills the context panel (Roles-10)", async () => {
      await session.step(98, "When user types \"BDD-RV-Role-{time}\" into gallery search", () => typeInto(page, session.text("BDD-RV-Role-{time}"), el("gallery search")));
      await session.step(99, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(100, "When user clicks on \"BDD-RV-Role-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-RV-Role-{time}\" link in gallery"))));
      await session.step(101, "Then the context panel should show \"BDD-RV-Role-{time}\"", () => contextPanelShows(page, session.text("BDD-RV-Role-{time}")));
      await session.step(102, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Actions\" accordion header in context panel"],["\"Assigned to\" accordion header in context panel"],["\"Favorites\" accordion header in context panel"],["\"Global Permissions\" accordion header in context panel"],["\"Permissions\" accordion header in context panel"]]));
      await session.step(108, "And MANAGE button in \"Assigned to\" section in context panel should be visible", () => shouldBe(page, el("MANAGE button in \"Assigned to\" section in context panel"), "visible"));
      await session.step(109, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(110, "Then no errors should have been logged", () => noErrors(page));
      await session.step(111, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
