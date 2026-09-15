/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/users-groups-roles/users-view.feature
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
import {clearField, clickOn, doubleClickOn, expand, followingShouldBe, pressKey, shouldBe, shouldHaveText, shouldNotBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, contextPanelOpen, contextPanelShows, dialogCloses, firstGalleryItemOther, galleryCountHigher, galleryCountLower, galleryCountNotLower, galleryMode, newUserOnServer, rememberFirstGalleryItem, rememberGalleryCount, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Users view", () => {
  const session = feature(test, "features/users-groups-roles/users-view.feature", import.meta.url);
  test("The Users view", {tag: ["@journey", "@users", "@realizes:views.users", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 16, page);
    await session.step(34, "Given user is logged in", () => loggedIn(page));
    await session.step(35, "And a new user \"opavlenko{time}v\" with email \"opavlenko+{time}v@datagrok.ai\" is on the server", () => newUserOnServer(page, session.text("opavlenko{time}v"), session.text("opavlenko+{time}v@datagrok.ai")));
    await session.step(36, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(37, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(38, "When user expands \"Platform\" tree node inside browse tree", () => expand(page, el("\"Platform\" tree node inside browse tree")));
    await session.step(39, "And user clicks on \"Platform > Users\" tree node inside browse tree", () => clickOn(page, el("\"Platform > Users\" tree node inside browse tree")));
    await run.scenario("The view opens from the Browse tree (Users-01)", async () => {
      await session.step(42, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(43, "And the page address should contain \"/users\"", () => urlShouldContain(page, "/users"));
      await session.step(44, "And gallery should be visible", () => shouldBe(page, el("gallery"), "visible"));
      await session.step(45, "And gallery counter should be visible", () => shouldBe(page, el("gallery counter"), "visible"));
      await session.step(46, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(47, "Then no errors should have been logged", () => noErrors(page));
      await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The toolbar carries its controls (Users-02)", async () => {
      await session.step(51, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["New button"],["gallery search"],["\"Switch to brief view\" icon inside gallery toolbar"],["\"Switch to card view\" icon inside gallery toolbar"],["\"Switch to grid view\" icon inside gallery toolbar"],["\"Sort list\" icon inside gallery toolbar"],["\"Toggle filters\" icon inside gallery toolbar"],["\"Refresh\" icon inside gallery toolbar"]]));
      await session.step(60, "And New button should be enabled", () => shouldBe(page, el("New button"), "enabled"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The view-mode icons switch the gallery's render mode (Users-11)", async () => {
      await session.step(65, "Then \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(66, "And the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(67, "When user clicks on \"Switch to card view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to card view\" icon inside gallery toolbar")));
      await session.step(68, "Then the gallery should be in card mode", () => galleryMode(page, "card"));
      await session.step(69, "And \"Switch to card view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(70, "And \"Switch to brief view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(71, "When user clicks on \"Switch to grid view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to grid view\" icon inside gallery toolbar")));
      await session.step(72, "Then the gallery should be in grid mode", () => galleryMode(page, "grid"));
      await session.step(73, "And \"Switch to grid view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to grid view\" icon inside gallery toolbar"), "selected"));
      await session.step(74, "And \"Switch to card view\" icon inside gallery toolbar should not be selected", () => shouldNotBe(page, el("\"Switch to card view\" icon inside gallery toolbar"), "selected"));
      await session.step(75, "When user clicks on \"Switch to brief view\" icon inside gallery toolbar", () => clickOn(page, el("\"Switch to brief view\" icon inside gallery toolbar")));
      await session.step(76, "Then the gallery should be in brief mode", () => galleryMode(page, "brief"));
      await session.step(77, "And \"Switch to brief view\" icon inside gallery toolbar should be selected", () => shouldBe(page, el("\"Switch to brief view\" icon inside gallery toolbar"), "selected"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Searching by login narrows the list, and clearing brings the rest back (Users-09)", async () => {
      await session.step(82, "When user types \"opavlenko{time}v\" into gallery search", () => typeInto(page, session.text("opavlenko{time}v"), el("gallery search")));
      await session.step(83, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(84, "And \"opavlenko{time}v\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"opavlenko{time}v\" link in gallery")), "visible"));
      await session.step(85, "And the page address should contain \"?q=opavlenko{time}v\"", () => urlShouldContain(page, session.text("?q=opavlenko{time}v")));
      await session.step(86, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(87, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(88, "Then the gallery counter should be higher than remembered", () => galleryCountHigher(page));
      await session.step(89, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(90, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(91, "Then no errors should have been logged", () => noErrors(page));
      await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("New offers a user, a service user and an invitation (Users-03)", async () => {
      await session.step(95, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(96, "Then the open menu should list \"User...\"", () => menuLists(page, "User..."));
      await session.step(97, "And the open menu should list \"Service User...\"", () => menuLists(page, "Service User..."));
      await session.step(98, "And the open menu should list \"Invite a Friend...\"", () => menuLists(page, "Invite a Friend..."));
      await session.step(99, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(100, "Then context menu should be hidden", () => shouldBe(page, el("context menu"), "hidden"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
      await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog asks for four fields, and Cancel closes it (Users-04)", async () => {
      await session.step(105, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(106, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(107, "Then \"Create new user\" dialog should be visible", () => shouldBe(page, el("\"Create new user\" dialog"), "visible"));
      await session.step(108, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["Email input in \"Create new user\" dialog"],["Login input in \"Create new user\" dialog"],["\"First Name\" input in \"Create new user\" dialog"],["\"Last Name\" input in \"Create new user\" dialog"]]));
      await session.step(113, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(114, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(115, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(116, "And no errors should have been logged", () => noErrors(page));
      await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The new user dialog refuses a bad email and a bad login (Users-06)", async () => {
      await session.step(120, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(121, "And user picks \"User...\" from the open menu", () => pickFromOpenMenu(page, "User..."));
      await session.step(122, "Then OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(123, "When user types \"not-an-email\" into Email input in \"Create new user\" dialog", () => typeInto(page, "not-an-email", el("Email input in \"Create new user\" dialog")));
      await session.step(124, "And user types \"Opavlenko+Bad\" into Login input in \"Create new user\" dialog", () => typeInto(page, "Opavlenko+Bad", el("Login input in \"Create new user\" dialog")));
      await session.step(125, "And user types \"Bad\" into \"First Name\" input in \"Create new user\" dialog", () => typeInto(page, "Bad", el("\"First Name\" input in \"Create new user\" dialog")));
      await session.step(126, "And user types \"Input\" into \"Last Name\" input in \"Create new user\" dialog", () => typeInto(page, "Input", el("\"Last Name\" input in \"Create new user\" dialog")));
      await session.step(127, "Then Email input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "invalid"));
      await session.step(128, "And Login input in \"Create new user\" dialog should be invalid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "invalid"));
      await session.step(129, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(130, "When user types \"opavlenko+bad{time}v@datagrok.ai\" into Email input in \"Create new user\" dialog", () => typeInto(page, session.text("opavlenko+bad{time}v@datagrok.ai"), el("Email input in \"Create new user\" dialog")));
      await session.step(131, "Then Email input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Email input in \"Create new user\" dialog"), "valid"));
      await session.step(132, "And OK button in \"Create new user\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "disabled"));
      await session.step(133, "When user types \"opavlenko-bad{time}v\" into Login input in \"Create new user\" dialog", () => typeInto(page, session.text("opavlenko-bad{time}v"), el("Login input in \"Create new user\" dialog")));
      await session.step(134, "Then Login input in \"Create new user\" dialog should be valid", () => shouldBe(page, el("Login input in \"Create new user\" dialog"), "valid"));
      await session.step(135, "And OK button in \"Create new user\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Create new user\" dialog"), "enabled"));
      await session.step(136, "When user clicks on CANCEL button in \"Create new user\" dialog", () => clickOn(page, el("CANCEL button in \"Create new user\" dialog")));
      await session.step(137, "Then the \"Create new user\" dialog should close", () => dialogCloses(page, "Create new user"));
      await session.step(138, "And no errors should have been logged", () => noErrors(page));
      await session.step(139, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Invite a Friend asks for an email (Users-08)", async () => {
      await session.step(142, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(143, "And user picks \"Invite a Friend...\" from the open menu", () => pickFromOpenMenu(page, "Invite a Friend..."));
      await session.step(144, "Then \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("\"Invite a Friend\" dialog"), "visible"));
      await session.step(145, "And Email input in \"Invite a Friend\" dialog should be visible", () => shouldBe(page, el("Email input in \"Invite a Friend\" dialog"), "visible"));
      await session.step(146, "When user clicks on CANCEL button in \"Invite a Friend\" dialog", () => clickOn(page, el("CANCEL button in \"Invite a Friend\" dialog")));
      await session.step(147, "Then the \"Invite a Friend\" dialog should close", () => dialogCloses(page, "Invite a Friend"));
      await session.step(148, "And no errors should have been logged", () => noErrors(page));
      await session.step(149, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A user's context menu (Users-14)", async () => {
      await session.step(152, "When user types \"opavlenko{time}v\" into gallery search", () => typeInto(page, session.text("opavlenko{time}v"), el("gallery search")));
      await session.step(153, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(154, "When user opens the context menu of \"opavlenko{time}v\" link in gallery", () => openContextMenu(page, el(session.text("\"opavlenko{time}v\" link in gallery"))));
      await session.step(155, "Then the open menu should list \"Details\"", () => menuLists(page, "Details"));
      await session.step(156, "And the open menu should list \"Chat\"", () => menuLists(page, "Chat"));
      await session.step(157, "And the open menu should list \"Disable...\"", () => menuLists(page, "Disable..."));
      await session.step(158, "And the open menu should list \"Groups...\"", () => menuLists(page, "Groups..."));
      await session.step(159, "And the open menu should list \"Roles...\"", () => menuLists(page, "Roles..."));
      await session.step(160, "And the open menu should list \"Copy > ID\"", () => menuLists(page, "Copy > ID"));
      await session.step(161, "And the open menu should list \"Add to favorites\"", () => menuLists(page, "Add to favorites"));
      await session.step(162, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(163, "And user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(164, "Then no errors should have been logged", () => noErrors(page));
      await session.step(165, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Selecting a user fills the context panel (Users-17)", async () => {
      await session.step(168, "When user types \"opavlenko{time}v\" into gallery search", () => typeInto(page, session.text("opavlenko{time}v"), el("gallery search")));
      await session.step(169, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(170, "When user clicks on \"opavlenko{time}v\" link in gallery", () => clickOn(page, el(session.text("\"opavlenko{time}v\" link in gallery"))));
      await session.step(171, "Then the context panel should show \"opavlenko{time}v\"", () => contextPanelShows(page, session.text("opavlenko{time}v")));
      await session.step(172, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Personal\" accordion header in context panel"],["\"Roles\" accordion header in context panel"],["\"Member of\" accordion header in context panel"],["\"Sticky meta\" accordion header in context panel"]]));
      await session.step(177, "And the following elements should be present:", () => followingShouldBe(page, "present", [["\"Projects\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Chats\" accordion header in context panel"],["\"Privileges\" accordion header in context panel"]]));
      await session.step(182, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(183, "Then no errors should have been logged", () => noErrors(page));
      await session.step(184, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Double-clicking a user opens the profile (Users-15)", async () => {
      await session.step(187, "When user types \"opavlenko{time}v\" into gallery search", () => typeInto(page, session.text("opavlenko{time}v"), el("gallery search")));
      await session.step(188, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(189, "When user double-clicks on \"opavlenko{time}v\" link in gallery", () => doubleClickOn(page, el(session.text("\"opavlenko{time}v\" link in gallery"))));
      await session.step(190, "Then the \"opavlenko{time}v\" view should be current", () => viewIsCurrent(page, session.text("opavlenko{time}v")));
      await session.step(191, "When user closes the current view", () => closeCurrentView(page));
      await session.step(192, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(193, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(194, "Then no errors should have been logged", () => noErrors(page));
      await session.step(195, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Details opens the same profile (Users-16)", async () => {
      await session.step(198, "When user types \"opavlenko{time}v\" into gallery search", () => typeInto(page, session.text("opavlenko{time}v"), el("gallery search")));
      await session.step(199, "Then the gallery counter should be lower than remembered", () => galleryCountLower(page));
      await session.step(200, "When user picks \"Details\" from the context menu of \"opavlenko{time}v\" link in gallery", () => pickFromContextMenu(page, "Details", el(session.text("\"opavlenko{time}v\" link in gallery"))));
      await session.step(201, "Then the \"opavlenko{time}v\" view should be current", () => viewIsCurrent(page, session.text("opavlenko{time}v")));
      await session.step(202, "When user closes the current view", () => closeCurrentView(page));
      await session.step(203, "Then the \"Users\" view should be current", () => viewIsCurrent(page, "Users"));
      await session.step(204, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(205, "Then no errors should have been logged", () => noErrors(page));
      await session.step(206, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A #tag search runs without an error (Users-10)", async () => {
      await session.step(209, "When user remembers the gallery counter", () => rememberGalleryCount(page));
      await session.step(210, "And user types \"#bddnotag{time}\" into gallery search", () => typeInto(page, session.text("#bddnotag{time}"), el("gallery search")));
      await session.step(211, "Then gallery counter should have text \"0\"", () => shouldHaveText(page, el("gallery counter"), "0"));
      await session.step(212, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(213, "Then the gallery counter should not be lower than remembered", () => galleryCountNotLower(page));
      await session.step(214, "And no errors should have been logged", () => noErrors(page));
      await session.step(215, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Sort list reorders the users, and Default reorders them back (Users-12)", async () => {
      await session.step(218, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(219, "And user clicks on \"Sort list\" icon inside gallery toolbar", () => clickOn(page, el("\"Sort list\" icon inside gallery toolbar")));
      await session.step(220, "Then the open menu should list \"Name\"", () => menuLists(page, "Name"));
      await session.step(221, "And the open menu should list \"Default\"", () => menuLists(page, "Default"));
      await session.step(222, "When user picks \"Name\" from the open menu", () => pickFromOpenMenu(page, "Name"));
      await session.step(223, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(224, "When user remembers the first item in gallery", () => rememberFirstGalleryItem(page));
      await session.step(225, "And user picks \"Default\" from the open menu", () => pickFromOpenMenu(page, "Default"));
      await session.step(226, "Then the first item in gallery should not be the remembered one", () => firstGalleryItemOther(page));
      await session.step(227, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(228, "Then no errors should have been logged", () => noErrors(page));
      await session.step(229, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Opening the filters for the first time logs no error (Users-13)", async () => {
      await session.step(238, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(239, "Then \"Joined\" filter card should be visible", () => shouldBe(page, el("\"Joined\" filter card"), "visible"));
      await session.step(240, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(241, "Then no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    await run.scenario("The filters show a card per user property (Users-13)", async () => {
      await session.step(244, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(245, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(246, "And \"Status\" filter card should be visible", () => shouldBe(page, el("\"Status\" filter card"), "visible"));
      await session.step(247, "And \"Joined\" filter card should be visible", () => shouldBe(page, el("\"Joined\" filter card"), "visible"));
      await session.step(248, "When user clicks on \"Toggle filters\" icon inside gallery toolbar", () => clickOn(page, el("\"Toggle filters\" icon inside gallery toolbar")));
      await session.step(249, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(250, "And no errors should have been logged", () => noErrors(page));
      await session.step(251, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
