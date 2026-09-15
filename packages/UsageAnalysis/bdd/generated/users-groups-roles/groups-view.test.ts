/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/groups-view.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.groups]
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
import {clearField, clickOn, expand, followingShouldBe, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, galleryCountHigher, galleryCountLower, galleryMode, groupOnServer, rememberGalleryCount, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Groups view", () => {
  const session = feature(test, "features/users-groups-roles/groups-view.feature", import.meta.url);
  test("The Groups view", {tag: ["@journey", "@groups", "@realizes:views.groups"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And a group named \"BDD-GV-Group-{time}\" is on the server", () => groupOnServer(page, session.text("BDD-GV-Group-{time}")));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(25, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(26, "And user clicks on \"Platform > Groups\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Groups\" tree node inside browse tree")));
    await run.scenario("The view opens from the Browse tree (Groups-01)", async () => {
      await session.step(29, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
      await session.step(30, "And the page address should contain \"/groups\"", () => urlShouldContain(page, "/groups"));
      await session.step(31, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(32, "And gallery counter should be visible", () => shouldBe(page, el("gallery counter"), "visible"));
      await session.step(33, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(34, "Then no errors should have been logged", () => noErrors(page));
      await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The toolbar carries its controls (Groups-02)", async () => {
      await session.step(38, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"New Group...\" button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]));
      await session.step(45, "And \"New Group...\" button should be enabled", () => shouldBe(page, el("\"New Group...\" button"), "enabled"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The view-mode icons switch the gallery's render mode (Groups-07)", async () => {
      await session.step(50, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(51, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(52, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(53, "And \"Switch to card view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(54, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(55, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(56, "And \"Switch to grid view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to grid view\" icon inside gallery toolbar"), "selected"));
      await session.step(57, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(58, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(59, "And \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Searching by name finds the one group, and clearing brings the rest back (Groups-06)", async () => {
      await session.step(64, "When user types \"BDD-GV-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GV-Group-{time}"), el("gallery search")));
      await session.step(65, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(66, "And \"BDD-GV-Group-{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDD-GV-Group-{time}\" link in gallery")), "visible"));
      await session.step(67, "And the page address should contain \"?q=BDD-GV-Group-{time}\"", () => urlShouldContain(page, session.text("?q=BDD-GV-Group-{time}")));
      await session.step(68, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(69, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(70, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(71, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(72, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(73, "Then no errors should have been logged", () => noErrors(page));
      await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A group's context menu (Groups-08)", async () => {
      await session.step(77, "When user types \"BDD-GV-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GV-Group-{time}"), el("gallery search")));
      await session.step(78, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(79, "When user opens the context menu of \"BDD-GV-Group-{time}\" link in gallery", () => openContextMenu(page, el(session.text("\"BDD-GV-Group-{time}\" link in gallery"))));
      await session.step(80, "Then the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(81, "And the open menu should list \"Request membership\"", () => menuLists(page, "Request membership"));
      await session.step(82, "And the open menu should list \"Chat\"", () => menuLists(page, "Chat"));
      await session.step(83, "And the open menu should list \"Delete\"", () => menuLists(page, "Delete"));
      await session.step(84, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(85, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(86, "Then no errors should have been logged", () => noErrors(page));
      await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting a group fills the context panel (Groups-10)", async () => {
      await session.step(90, "When user types \"BDD-GV-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GV-Group-{time}"), el("gallery search")));
      await session.step(91, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(92, "When user clicks on \"BDD-GV-Group-{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDD-GV-Group-{time}\" link in gallery"))));
      await session.step(93, "Then the context panel should show \"BDD-GV-Group-{time}\"", () => contextPanelShows(page, session.text("BDD-GV-Group-{time}")));
      await session.step(94, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Actions\" accordion header in context panel"],["\"Members\" accordion header in context panel"],["\"Favorites\" accordion header in context panel"],["\"Global Permissions\" accordion header in context panel"],["\"Permissions\" accordion header in context panel"],["\"Sticky meta\" accordion header in context panel"]]));
      await session.step(101, "And MANAGE button in \"Members\" section in context panel should be visible", () => shouldBe(page, el("MANAGE button in \"Members\" section in context panel"), "visible"));
      await session.step(102, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(103, "Then no errors should have been logged", () => noErrors(page));
      await session.step(104, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Chat opens the group's own chat (Groups-17)", async () => {
      await session.step(107, "When user types \"BDD-GV-Group-{time}\" into gallery search", () => typeInto(page, session.text("BDD-GV-Group-{time}"), el("gallery search")));
      await session.step(108, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(109, "When user picks \"Chat\" from the context menu of \"BDD-GV-Group-{time}\" link in gallery", () => pickFromContextMenu(page, "Chat", el(session.text("\"BDD-GV-Group-{time}\" link in gallery"))));
      await session.step(110, "Then the \"Chats\" view should be current", () => viewIsCurrent(page, "Chats"));
      await session.step(111, "And chat header should contain text \"BDD-GV-Group-{time}\"", () => shouldContainText(page, el("chat header"), session.text("BDD-GV-Group-{time}")));
      await session.step(112, "When user closes the current view", () => closeCurrentView(page));
      await session.step(113, "Then the \"Groups\" view should be current", () => viewIsCurrent(page, "Groups"));
      await session.step(114, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(115, "Then no errors should have been logged", () => noErrors(page));
      await session.step(116, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
