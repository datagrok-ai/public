/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-view.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.users]
--- */
import {test} from '@playwright/test';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, expand, followingShouldBe, pressKey, shouldBe, shouldHaveText, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, dialogCloses, firstGalleryItemOther, galleryCountHigher, galleryCountLower, galleryCountNotLower, galleryMode, rememberFirstGalleryItem, rememberGalleryCount, urlShouldContain, urlShouldNotContain, userOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Users view", () => {
  const session = feature(test, "features/users-groups-roles/users-view.feature", import.meta.url);
  test("The Users view", {tag: ["@journey", "@users", "@realizes:views.users"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 15, page);
    await session.step(33, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "And a user \"bddviewed\" is on the server", () => userOnServer(page, "bddviewed"));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(37, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(38, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await run.scenario("The view opens from the Browse tree (Users-01)", async () => {
      await session.step(41, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(42, "And the page address should contain \"/users\"", () => urlShouldContain(page, "/users"));
      await session.step(43, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(44, "And gallery counter should be visible", () => shouldBe(page, el("gallery counter"), "visible"));
      await session.step(45, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(46, "Then no errors should have been logged", () => noErrors(page));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The toolbar carries its controls (Users-02)", async () => {
      await session.step(50, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["New button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Sort list\" icon inside gallery toolbar"],["\"Toggle filters\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]), [["New button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Sort list\" icon inside gallery toolbar"],["\"Toggle filters\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]);
      await session.step(59, "And New button should be enabled", () => shouldBe(page, el("New button"), "enabled"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
      await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The view-mode icons switch the gallery's render mode (Users-11)", async () => {
      await session.step(64, "Then \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(65, "And the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(66, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(67, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(68, "And \"Switch to card view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(69, "And \"Switch to brief view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(70, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(71, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(72, "And \"Switch to grid view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to grid view\" icon inside gallery toolbar"), "selected"));
      await session.step(73, "And \"Switch to card view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(74, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(75, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(76, "And \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Searching by login narrows the list, and clearing brings the rest back (Users-09)", async () => {
      await session.step(81, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(82, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(83, "And \"bddviewed\" link in gallery should be visible", () => shouldBe(page, el("\"bddviewed\" link in gallery"), "visible"));
      await session.step(84, "And the page address should contain \"?q=bddviewed\"", () => urlShouldContain(page, "?q=bddviewed"));
      await session.step(85, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(86, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(87, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(88, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(89, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(90, "Then no errors should have been logged", () => noErrors(page));
      await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("New offers a user, a service user and an invitation (Users-03)", async () => {
      await session.step(94, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(95, "Then the open menu should list \"User...\"", () => menuLists(page, "User..."));
      await session.step(96, "And the open menu should list \"Service User...\"", () => menuLists(page, "Service User..."));
      await session.step(97, "And the open menu should list \"Invite a Friend...\"", () => menuLists(page, "Invite a Friend..."));
      await session.step(98, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(99, "Then context menu should be hidden", () => shouldBe(page, el("context menu"), "hidden"));
      await session.step(100, "And no errors should have been logged", () => noErrors(page));
      await session.step(101, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog asks for four fields, and Cancel closes it (Users-04)", async () => {
      await session.step(104, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(105, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(106, "Then \"Create new user\" dialog should be visible", () => shouldBe(page, el("\"Create new user\" dialog"), "visible"));
      await session.step(107, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Email input in \"Create new user\" dialog"],["Login input in \"Create new user\" dialog"],["\"First Name\" input in \"Create new user\" dialog"],["\"Last Name\" input in \"Create new user\" dialog"]]), [["Email input in \"Create new user\" dialog"],["Login input in \"Create new user\" dialog"],["\"First Name\" input in \"Create new user\" dialog"],["\"Last Name\" input in \"Create new user\" dialog"]]);
      await session.step(112, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(113, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(114, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
      await session.step(116, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog refuses a bad email and a bad login (Users-06)", async () => {
      await session.step(119, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(120, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(121, "Then OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(122, "When user types \"not-an-email\" into Email input in \"Create new user\" dialog", () => typeInto(page, "not-an-email", el("Email input in \"Create new user\" dialog")));
      await session.step(123, "And user types \"Bdd+Bad\" into Login input in \"Create new user\" dialog", () => typeInto(page, "Bdd+Bad", el("Login input in \"Create new user\" dialog")));
      await session.step(124, "And user types \"Bad\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "Bad", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(125, "And user types \"Input\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, "Input", el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(126, "Then Email input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "invalid"));
      await session.step(127, "And Login input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "invalid"));
      await session.step(128, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(129, "When user types \"bdd-bad{time}@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd-bad{time}@datagrok.ai"), el("Email input in \"Create new user\" dialog")));
      await session.step(130, "Then Email input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "valid"));
      await session.step(131, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(132, "When user types \"bdd-bad{time}\" into Login input in \"Create new user\" dialog", () => typeInto(page, session.text("bdd-bad{time}"), el("Login input in \"Create new user\" dialog")));
      await session.step(133, "Then Login input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "valid"));
      await session.step(134, "And OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(135, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(136, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
      await session.step(138, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Invite a Friend asks for an email (Users-08)", async () => {
      await session.step(141, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(142, "And user picks \"Invite a Friend...\" from the open menu", () => pickFromOpenMenu(page, "Invite a Friend..."));
      await session.step(143, "Then \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("\"Invite a Friend\" dialog"), "visible"));
      await session.step(144, "And Email input in \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("Email input in \"Invite a Friend\" dialog"), "visible"));
      await session.step(145, "When user clicks on CANCEL button in \"Invite a Friend\" dialog", () => clickOn(page, el("CANCEL button in \"Invite a Friend\" dialog")));
      await session.step(146, "Then the \"Invite a Friend\" dialog should close", () => dialogCloses(page, "Invite a Friend"));
      await session.step(147, "And no errors should have been logged", () => noErrors(page));
      await session.step(148, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's context menu (Users-14)", async () => {
      await session.step(151, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(152, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(153, "When user opens the context menu of \"bddviewed\" link in gallery", () => openContextMenu(page, el("\"bddviewed\" link in gallery")));
      await session.step(154, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
      await session.step(155, "And the open menu should list \"Chat\"", () => menuLists(page, "Chat"));
      await session.step(156, "And the open menu should list \"Disable...\"", () => menuLists(page, "Disable..."));
      await session.step(157, "And the open menu should list \"Groups...\"", () => menuLists(page, "Groups..."));
      await session.step(158, "And the open menu should list \"Roles...\"", () => menuLists(page, "Roles..."));
      await session.step(159, "And the open menu should list \"Copy > ID\"", () => menuLists(page, "Copy > ID"));
      await session.step(160, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
      await session.step(161, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(162, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(163, "Then no errors should have been logged", () => noErrors(page));
      await session.step(164, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting a user fills the context panel (Users-17)", async () => {
      await session.step(167, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(168, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(169, "When user clicks on \"bddviewed\" link in gallery", () => clickOn(page, el("\"bddviewed\" link in gallery")));
      await session.step(170, "Then the context panel should show \"bddviewed\"", () => contextPanelShows(page, "bddviewed"));
      await session.step(171, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Personal\" accordion header in context panel"],["\"Roles\" accordion header in context panel"],["\"Member of\" accordion header in context panel"]]), [["\"Personal\" accordion header in context panel"],["\"Roles\" accordion header in context panel"],["\"Member of\" accordion header in context panel"]]);
      await session.step(175, "And the following elements should be present:", () => followingShouldBe(page, "present", [["\"Projects\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Chats\" accordion header in context panel"],["\"Privileges\" accordion header in context panel"]]), [["\"Projects\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Chats\" accordion header in context panel"],["\"Privileges\" accordion header in context panel"]]);
      await session.step(180, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(181, "Then no errors should have been logged", () => noErrors(page));
      await session.step(182, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Double-clicking a user opens the profile (Users-15)", async () => {
      await session.step(185, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(186, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(187, "When user double-clicks on \"bddviewed\" link in gallery", () => doubleClickOn(page, el("\"bddviewed\" link in gallery")));
      await session.step(188, "Then the \"bddviewed\" view should be current", () => viewIsCurrent(page, "bddviewed"));
      await session.step(189, "When user closes the current view", () => closeCurrentView(page));
      await session.step(190, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(191, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(192, "Then no errors should have been logged", () => noErrors(page));
      await session.step(193, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Details opens the same profile (Users-16)", async () => {
      await session.step(196, "When user types \"bddviewed\" into gallery search", () => typeInto(page, "bddviewed", el("gallery search")));
      await session.step(197, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(198, "When user picks \"Details\" from the context menu of \"bddviewed\" link in gallery", () => pickFromContextMenu(page, "Details", el("\"bddviewed\" link in gallery")));
      await session.step(199, "Then the \"bddviewed\" view should be current", () => viewIsCurrent(page, "bddviewed"));
      await session.step(200, "When user closes the current view", () => closeCurrentView(page));
      await session.step(201, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(202, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(203, "Then no errors should have been logged", () => noErrors(page));
      await session.step(204, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A #tag search runs without an error (Users-10)", async () => {
      await session.step(207, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(208, "And user types \"#bddnotag{time}\" into gallery search", () => typeInto(page, session.text("#bddnotag{time}"), el("gallery search")));
      await session.step(209, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(210, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(211, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(212, "And no errors should have been logged", () => noErrors(page));
      await session.step(213, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Sort list reorders the users, and Default reorders them back (Users-12)", async () => {
      await session.step(216, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(217, "And user clicks on \"Sort list\" icon inside gallery toolbar", () => clickOn(page, el("\"Sort list\" icon inside gallery toolbar")));
      await session.step(218, "Then the open menu should list \"Name\"", () => menuLists(page, "Name"));
      await session.step(219, "And the open menu should list \"Default\"", () => menuLists(page, "Default"));
      await session.step(220, "When user picks \"Name\" from the open menu", () => pickFromOpenMenu(page, "Name"));
      await session.step(221, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(222, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(223, "And user picks \"Default\" from the open menu", () => pickFromOpenMenu(page, "Default"));
      await session.step(224, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(225, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(226, "Then no errors should have been logged", () => noErrors(page));
      await session.step(227, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The filters show a card per user property (Users-13)", async () => {
      await session.step(230, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(231, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(232, "And \"Status\" filter card should be visible", () => shouldBe(page, el("\"Status\" filter card"), "visible"));
      await session.step(233, "And \"Joined\" filter card should be visible", () => shouldBe(page, el("\"Joined\" filter card"), "visible"));
      await session.step(234, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(235, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(236, "And no errors should have been logged", () => noErrors(page));
      await session.step(237, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
