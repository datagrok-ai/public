/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-view.feature
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
import {clearField, clickOn, doubleClickOn, enterInto, expand, followingShouldBe, pressKey, shouldBe, shouldContainText, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, dialogCloses, firstGalleryItemOther, galleryCountHigher, galleryCountLower, galleryMode, rememberFirstGalleryItem, rememberGalleryCount, urlShouldContain, urlShouldNotContain, userOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuDoesNotList, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Users view", () => {
  const session = feature(test, "features/users-groups-roles/users-view.feature", import.meta.url);
  test("The Users view", {tag: ["@journey", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(39, "Given user is logged in", () => loggedIn(page));
    await session.step(40, "And a user \"bddviewed\" is on the server", () => userOnServer(page, "bddviewed"));
    await session.step(41, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(42, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(43, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(44, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await run.scenario("The view opens from the Browse tree (Users-01)", async () => {
      await session.step(47, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(48, "And the page address should contain \"/users\"", () => urlShouldContain(page, "/users"));
      await session.step(49, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(50, "And gallery counter should be visible", () => shouldBe(page, el("gallery counter"), "visible"));
      await session.step(51, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
      await session.step(53, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The toolbar carries its controls (Users-02)", async () => {
      await session.step(56, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["New button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Sort list\" icon inside gallery toolbar"],["\"Toggle filters\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]), [["New button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Sort list\" icon inside gallery toolbar"],["\"Toggle filters\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]);
      await session.step(65, "And New button should be enabled", () => shouldBe(page, el("New button"), "enabled"));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The view-mode icons switch the gallery's render mode (Users-11)", async () => {
      await session.step(70, "Then \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(71, "And the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(72, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(73, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(74, "And \"Switch to card view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(75, "And \"Switch to brief view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(76, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(77, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(78, "And \"Switch to grid view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to grid view\" icon inside gallery toolbar"), "selected"));
      await session.step(79, "And \"Switch to card view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(80, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(81, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(82, "And \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
      await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Searching by login narrows the list, and clearing brings the rest back (Users-09)", async () => {
      await session.step(87, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(88, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(89, "And \"bddviewed\" link in gallery should be visible", () => shouldBe(page, el("\"bddviewed\" link in gallery"), "visible"));
      await session.step(90, "And the page address should contain \"?q=bddviewed\"", () => urlShouldContain(page, "?q=bddviewed"));
      await session.step(91, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(92, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(93, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(94, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(95, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(96, "Then no errors should have been logged", () => noErrors(page));
      await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("New offers a user, a service user and an invitation (Users-03)", async () => {
      await session.step(100, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(101, "Then the open menu should list \"User...\"", () => menuLists(page, "User..."));
      await session.step(102, "And the open menu should list \"Service User...\"", () => menuLists(page, "Service User..."));
      await session.step(103, "And the open menu should list \"Invite a Friend...\"", () => menuLists(page, "Invite a Friend..."));
      await session.step(104, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(105, "Then context menu should be hidden", () => shouldBe(page, el("context menu"), "hidden"));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
      await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog asks for four fields, and Cancel closes it (Users-04)", async () => {
      await session.step(110, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(111, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(112, "Then \"Create new user\" dialog should be visible", () => shouldBe(page, el("\"Create new user\" dialog"), "visible"));
      await session.step(113, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Email input in \"Create new user\" dialog"],["Login input in \"Create new user\" dialog"],["\"First Name\" input in \"Create new user\" dialog"],["\"Last Name\" input in \"Create new user\" dialog"]]), [["Email input in \"Create new user\" dialog"],["Login input in \"Create new user\" dialog"],["\"First Name\" input in \"Create new user\" dialog"],["\"Last Name\" input in \"Create new user\" dialog"]]);
      await session.step(118, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(119, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(120, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
      await session.step(122, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog refuses a bad email and a bad login (Users-06)", async () => {
      await session.step(125, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(126, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(127, "Then OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(128, "When user types \"not-an-email\" into Email input in \"Create new user\" dialog", () => typeInto(page, "not-an-email", el("Email input in \"Create new user\" dialog")));
      await session.step(129, "And user types \"Bdd+Bad\" into Login input in \"Create new user\" dialog", () => typeInto(page, "Bdd+Bad", el("Login input in \"Create new user\" dialog")));
      await session.step(130, "And user types \"Bad\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "Bad", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(131, "And user types \"Input\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, "Input", el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(132, "Then Email input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "invalid"));
      await session.step(133, "And Login input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "invalid"));
      await session.step(134, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(135, "When user types \"bdd-bad{time}@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd-bad{time}@datagrok.ai"), el("Email input in \"Create new user\" dialog")));
      await session.step(136, "Then Email input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "valid"));
      await session.step(137, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(138, "When user types \"bdd-bad{time}\" into Login input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd-bad{time}"), el("Login input in \"Create new user\" dialog")));
      await session.step(139, "Then Login input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "valid"));
      await session.step(140, "And OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(141, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(142, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
      await session.step(144, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Invite a Friend asks for an email (Users-08)", async () => {
      await session.step(147, "Given a user \"bddprofile\" is on the server", () => userOnServer(page, "bddprofile"));
      await session.step(148, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(149, "And user picks \"Invite a Friend...\" from the open menu", () => pickFromOpenMenu(page, "Invite a Friend..."));
      await session.step(150, "Then \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("\"Invite a Friend\" dialog"), "visible"));
      await session.step(151, "And Email input in \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("Email input in \"Invite a Friend\" dialog"), "visible"));
      await session.step(152, "When user enters \"bddprofile@datagrok.ai\" into Email input in \"Invite a Friend\" dialog", () => enterInto(page, "bddprofile@datagrok.ai", el("Email input in \"Invite a Friend\" dialog")));
      await session.step(153, "Then \"Invite a Friend\" dialog should contain text \"User is already registered\"", () => shouldContainText(page, el("\"Invite a Friend\" dialog"), "User is already registered"));
      await session.step(154, "When user clicks on CANCEL button in \"Invite a Friend\" dialog", () => clickOn(page, el("CANCEL button in \"Invite a Friend\" dialog")));
      await session.step(155, "Then the \"Invite a Friend\" dialog should close", () => dialogCloses(page, "Invite a Friend"));
      await session.step(156, "And no errors should have been logged", () => noErrors(page));
      await session.step(157, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's context menu (Users-14)", async () => {
      await session.step(160, "When user types \"{login}\" into gallery search", () => typeInto(page, session.text("{login}"), el("gallery search")));
      await session.step(161, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(162, "When user opens the context menu of \"{login}\" link in gallery", () => openContextMenu(page, el(session.text("\"{login}\" link in gallery"))));
      await session.step(163, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
      await session.step(164, "And the open menu should not list \"Disable...\"", () => menuDoesNotList(page, "Disable..."));
      await session.step(165, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(166, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(167, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(168, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(169, "When user opens the context menu of \"bddviewed\" link in gallery", () => openContextMenu(page, el("\"bddviewed\" link in gallery")));
      await session.step(170, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
      await session.step(171, "And the open menu should list \"Chat\"", () => menuLists(page, "Chat"));
      await session.step(172, "And the open menu should list \"Disable...\"", () => menuLists(page, "Disable..."));
      await session.step(173, "And the open menu should list \"Groups...\"", () => menuLists(page, "Groups..."));
      await session.step(174, "And the open menu should list \"Roles...\"", () => menuLists(page, "Roles..."));
      await session.step(175, "And the open menu should list \"Copy > ID\"", () => menuLists(page, "Copy > ID"));
      await session.step(176, "And the open menu should list \"Add To Favorites\"", () => menuLists(page, "Add To Favorites"));
      await session.step(177, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(178, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(179, "Then no errors should have been logged", () => noErrors(page));
      await session.step(180, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting a user fills the context panel (Users-17)", async () => {
      await session.step(183, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(184, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(185, "When user clicks on \"bddviewed\" link in gallery", () => clickOn(page, el("\"bddviewed\" link in gallery")));
      await session.step(186, "Then the context panel should show \"bddviewed\"", () => contextPanelShows(page, "bddviewed"));
      await session.step(187, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Personal\" accordion header in context panel"],["\"Roles\" accordion header in context panel"],["\"Member of\" accordion header in context panel"]]), [["\"Personal\" accordion header in context panel"],["\"Roles\" accordion header in context panel"],["\"Member of\" accordion header in context panel"]]);
      await session.step(191, "And the following elements should be present:", () => followingShouldBe(page, "present", [["\"Projects\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Chats\" accordion header in context panel"],["\"Privileges\" accordion header in context panel"]]), [["\"Projects\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Chats\" accordion header in context panel"],["\"Privileges\" accordion header in context panel"]]);
      await session.step(196, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(197, "Then no errors should have been logged", () => noErrors(page));
      await session.step(198, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Double-clicking a user opens the profile (Users-15)", async () => {
      await session.step(201, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(202, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(203, "When user double-clicks on \"bddviewed\" link in gallery", () => doubleClickOn(page, el("\"bddviewed\" link in gallery")));
      await session.step(204, "Then the \"bddviewed\" view should be current", () => viewIsCurrent(page, "bddviewed"));
      await session.step(205, "When user closes the current view", () => closeCurrentView(page));
      await session.step(206, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(207, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(208, "Then no errors should have been logged", () => noErrors(page));
      await session.step(209, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Details opens the same profile (Users-16)", async () => {
      await session.step(212, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(213, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(214, "When user picks \"Details\" from the context menu of \"bddviewed\" link in gallery", () => pickFromContextMenu(page, "Details", el("\"bddviewed\" link in gallery")));
      await session.step(215, "Then the \"bddviewed\" view should be current", () => viewIsCurrent(page, "bddviewed"));
      await session.step(216, "When user closes the current view", () => closeCurrentView(page));
      await session.step(217, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(218, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(219, "Then no errors should have been logged", () => noErrors(page));
      await session.step(220, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Sort list reorders the users, and Default reorders them back (Users-12)", async () => {
      await session.step(223, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(224, "And user clicks on \"Sort list\" icon inside gallery toolbar", () => clickOn(page, el("\"Sort list\" icon inside gallery toolbar")));
      await session.step(225, "Then the open menu should list \"Name\"", () => menuLists(page, "Name"));
      await session.step(226, "And the open menu should list \"Default\"", () => menuLists(page, "Default"));
      await session.step(227, "When user picks \"Name\" from the open menu", () => pickFromOpenMenu(page, "Name"));
      await session.step(228, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(229, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(230, "And user picks \"Default\" from the open menu", () => pickFromOpenMenu(page, "Default"));
      await session.step(231, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(232, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(233, "Then no errors should have been logged", () => noErrors(page));
      await session.step(234, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The filters show a card per user property (Users-13)", async () => {
      await session.step(237, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(238, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(239, "And \"Status\" filter card should be visible", () => shouldBe(page, el("\"Status\" filter card"), "visible"));
      await session.step(240, "And \"Joined\" filter card should be visible", () => shouldBe(page, el("\"Joined\" filter card"), "visible"));
      await session.step(241, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(242, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(243, "And no errors should have been logged", () => noErrors(page));
      await session.step(244, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
