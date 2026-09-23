/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/home/home-widgets.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.view.welcome, powerpack.dashboard.spotlight, powerpack.dashboard.community]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {adminAccountOnServer, containsOneOf, everyWidgetHasContent, homeWidgetsAre, linksPointTo, placeholderStarts, reloadPage, signedInAs, tabShowing, tipOfTheDay, viewOpen, widgetSettingsRestored, widgetStoredAs} from '../../bindings/home.js';
import {searchFinished} from '../../bindings/search.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, followingShouldBe, hoverOver, shouldBe, shouldContainText, shouldHaveText, shouldNotBe, typeInto, uncheck, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, openDataset, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The widgets of the Home page", () => {
  const session = feature(test, "features/home/home-widgets.feature", import.meta.url);
  test("The widgets of the Home page", {tag: ["@journey", "@serial", "@realizes:powerpack.view.welcome", "@realizes:powerpack.dashboard.spotlight", "@realizes:powerpack.dashboard.community", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(40, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "And an administrator account \"bddhomeadmin\" is on the server", () => adminAccountOnServer(page, "bddhomeadmin"));
    await session.step(42, "And user is signed in as \"bddhomeadmin\" on this page", () => signedInAs(page, "bddhomeadmin"));
    await session.step(43, "And the widget settings of the Home page come back when the feature ends", () => widgetSettingsRestored(page));
    await run.scenario("The Home page shows the search box and four loaded widgets after a reload", async () => {
      await session.step(46, "When user reloads the page", () => reloadPage(page));
      await session.step(47, "Then home search should be visible", () => shouldBe(page, el("home search"), "visible"));
      await session.step(48, "And the placeholder of home search should start with \"Search everywhere\"", () => placeholderStarts(page, el("home search"), "Search everywhere"));
      await session.step(49, "And the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(50, "And every widget of the Home page should show content", () => everyWidgetHasContent(page));
      await session.step(51, "And title of Community home widget should have text \"Community\"", () => shouldHaveText(page, el("title of Community home widget"), "Community"));
      await session.step(52, "And title of Spotlight home widget should be hidden", () => shouldBe(page, el("title of Spotlight home widget"), "hidden"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The close icon shows on hover and removes the widget", async () => {
      await session.step(57, "Then close icon of Community home widget should be hidden", () => shouldBe(page, el("close icon of Community home widget"), "hidden"));
      await session.step(58, "When user hovers over Community home widget", () => hoverOver(page, el("Community home widget")));
      await session.step(59, "Then close icon of Community home widget should be visible", () => shouldBe(page, el("close icon of Community home widget"), "visible"));
      await session.step(60, "When user clicks on close icon of Community home widget", () => clickOn(page, el("close icon of Community home widget")));
      await session.step(61, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(62, "And the Home page should show the widgets \"Spotlight, Reports, Usage\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage"));
      await session.step(63, "And the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(64, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(65, "And user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(66, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(67, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A closed widget stays closed after a reload", async () => {
      await session.step(71, "When user hovers over Community home widget", () => hoverOver(page, el("Community home widget")));
      await session.step(72, "And user clicks on close icon of Community home widget", () => clickOn(page, el("close icon of Community home widget")));
      await session.step(73, "Then the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(74, "When user reloads the page", () => reloadPage(page));
      await session.step(75, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(76, "And the Home page should show the widgets \"Spotlight, Reports, Usage\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage"));
      await session.step(77, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(78, "Then \"Community\" input in context panel should not be checked", () => shouldNotBe(page, el("\"Community\" input in context panel"), "checked"));
      await session.step(79, "When user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(80, "Then the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(81, "When user reloads the page", () => reloadPage(page));
      await session.step(82, "Then the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Customize form hides and shows a widget", async () => {
      await session.step(86, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(87, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Spotlight\" input in context panel"],["\"Community\" input in context panel"],["\"Usage\" input in context panel"],["\"Reports\" input in context panel"]]), [["\"Spotlight\" input in context panel"],["\"Community\" input in context panel"],["\"Usage\" input in context panel"],["\"Reports\" input in context panel"]]);
      await session.step(92, "And \"Community\" input in context panel should be checked", () => shouldBe(page, el("\"Community\" input in context panel"), "checked"));
      await session.step(93, "When user unchecks \"Community\" input in context panel", () => uncheck(page, el("\"Community\" input in context panel")));
      await session.step(94, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(95, "And the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(96, "When user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(97, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(98, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A widget hidden in the Customize form stays hidden after a reload", async () => {
      await session.step(102, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(103, "And user unchecks \"Community\" input in context panel", () => uncheck(page, el("\"Community\" input in context panel")));
      await session.step(104, "Then the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(105, "When user reloads the page", () => reloadPage(page));
      await session.step(106, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(107, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(108, "Then \"Community\" input in context panel should not be checked", () => shouldNotBe(page, el("\"Community\" input in context panel"), "checked"));
      await session.step(109, "When user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(110, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(111, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A search replaces the widgets, and clearing it brings them back", async () => {
      await session.step(115, "When user types \"aspirin\" into home search", () => typeInto(page, "aspirin", el("home search")));
      await session.step(116, "Then home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
      await session.step(117, "And home search results should be visible", () => shouldBe(page, el("home search results"), "visible"));
      await session.step(118, "And the page address should contain \"search?q=aspirin\"", () => urlShouldContain(page, "search?q=aspirin"));
      await session.step(119, "And the search should have finished", () => searchFinished(page));
      await session.step(120, "When user clears home search", () => clearField(page, el("home search")));
      await session.step(121, "Then home widgets panel should be visible", () => shouldBe(page, el("home widgets panel"), "visible"));
      await session.step(122, "And home search results should be hidden", () => shouldBe(page, el("home search results"), "hidden"));
      await session.step(123, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Home icon brings back the same widgets, before and after a reload", async () => {
      await session.step(127, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(128, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(129, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(130, "When user clicks on \"Home\" icon in browse toolbar", () => clickOn(page, el("\"Home\" icon in browse toolbar")));
      await session.step(131, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(132, "And the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(133, "When user reloads the page", () => reloadPage(page));
      await session.step(134, "Then the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(135, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Spotlight has six tabs, and each shows its own page", async () => {
      await session.step(138, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Workspace tab in Spotlight home widget"],["Spotlight tab in Spotlight home widget"],["Favorites tab in Spotlight home widget"],["Notifications tab in Spotlight home widget"],["\"My Activity\" tab in Spotlight home widget"],["Learn tab in Spotlight home widget"]]), [["Workspace tab in Spotlight home widget"],["Spotlight tab in Spotlight home widget"],["Favorites tab in Spotlight home widget"],["Notifications tab in Spotlight home widget"],["\"My Activity\" tab in Spotlight home widget"],["Learn tab in Spotlight home widget"]]);
      await session.step(145, "When user clicks on Spotlight tab in Spotlight home widget", () => clickOn(page, el("Spotlight tab in Spotlight home widget")));
      await session.step(146, "Then the \"Spotlight\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Spotlight", el("Spotlight home widget")));
      await session.step(147, "And spotlight page of Spotlight home widget should contain one of the texts \"Recent | Interactive Tutorials\"", () => containsOneOf(page, el("spotlight page of Spotlight home widget"), "Recent | Interactive Tutorials"));
      await session.step(148, "When user clicks on Favorites tab in Spotlight home widget", () => clickOn(page, el("Favorites tab in Spotlight home widget")));
      await session.step(149, "Then the \"Favorites\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Favorites", el("Spotlight home widget")));
      await session.step(150, "When user clicks on Notifications tab in Spotlight home widget", () => clickOn(page, el("Notifications tab in Spotlight home widget")));
      await session.step(151, "Then the \"Notifications\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Notifications", el("Spotlight home widget")));
      await session.step(152, "When user clicks on \"My Activity\" tab in Spotlight home widget", () => clickOn(page, el("\"My Activity\" tab in Spotlight home widget")));
      await session.step(153, "Then the \"My Activity\" tab of Spotlight home widget should be showing", () => tabShowing(page, "My Activity", el("Spotlight home widget")));
      await session.step(154, "When user clicks on Learn tab in Spotlight home widget", () => clickOn(page, el("Learn tab in Spotlight home widget")));
      await session.step(155, "Then the \"Learn\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Learn", el("Spotlight home widget")));
      await session.step(156, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["VIDEO tab in Spotlight home widget"],["WIKI tab in Spotlight home widget"],["DEMO tab in Spotlight home widget"],["TUTORIALS tab in Spotlight home widget"]]), [["VIDEO tab in Spotlight home widget"],["WIKI tab in Spotlight home widget"],["DEMO tab in Spotlight home widget"],["TUTORIALS tab in Spotlight home widget"]]);
      await session.step(161, "And Spotlight home widget should contain text \"Cheminformatics\"", () => shouldContainText(page, el("Spotlight home widget"), "Cheminformatics"));
      await session.step(162, "When user clicks on WIKI tab in Spotlight home widget", () => clickOn(page, el("WIKI tab in Spotlight home widget")));
      await session.step(163, "Then the \"WIKI\" tab of Spotlight home widget should be showing", () => tabShowing(page, "WIKI", el("Spotlight home widget")));
      await session.step(164, "When user clicks on DEMO tab in Spotlight home widget", () => clickOn(page, el("DEMO tab in Spotlight home widget")));
      await session.step(165, "Then the \"DEMO\" tab of Spotlight home widget should be showing", () => tabShowing(page, "DEMO", el("Spotlight home widget")));
      await session.step(166, "When user clicks on TUTORIALS tab in Spotlight home widget", () => clickOn(page, el("TUTORIALS tab in Spotlight home widget")));
      await session.step(167, "Then the \"TUTORIALS\" tab of Spotlight home widget should be showing", () => tabShowing(page, "TUTORIALS", el("Spotlight home widget")));
      await session.step(168, "When user clicks on VIDEO tab in Spotlight home widget", () => clickOn(page, el("VIDEO tab in Spotlight home widget")));
      await session.step(169, "Then the \"VIDEO\" tab of Spotlight home widget should be showing", () => tabShowing(page, "VIDEO", el("Spotlight home widget")));
      await session.step(170, "When user clicks on Workspace tab in Spotlight home widget", () => clickOn(page, el("Workspace tab in Spotlight home widget")));
      await session.step(171, "Then the \"Workspace\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Workspace", el("Spotlight home widget")));
      await session.step(172, "And no errors should have been logged", () => noErrors(page));
      await session.step(173, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The tip of the day at the bottom of Spotlight opens what it names", async () => {
      await session.step(176, "Then the tip of the day of Spotlight home widget should open what it names", () => tipOfTheDay(page, el("Spotlight home widget")));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Community lists links to the community site", async () => {
      await session.step(180, "Then every link of Community home widget should point to \"https://community.datagrok.ai/\"", () => linksPointTo(page, el("Community home widget"), "https://community.datagrok.ai/"));
      await session.step(181, "And Community home widget should contain text \"Platform Releases\"", () => shouldContainText(page, el("Community home widget"), "Platform Releases"));
    });
    await run.scenario("Usage shows its user and error charts and the state of the services", async () => {
      await session.step(184, "Then Usage home widget should contain text \"Users\"", () => shouldContainText(page, el("Usage home widget"), "Users"));
      await session.step(185, "And Usage home widget should contain text \"Errors\"", () => shouldContainText(page, el("Usage home widget"), "Errors"));
      await session.step(186, "And Usage home widget should contain text \"System\"", () => shouldContainText(page, el("Usage home widget"), "System"));
      await session.step(187, "And Usage home widget should contain text \"Jupyter\"", () => shouldContainText(page, el("Usage home widget"), "Jupyter"));
      await session.step(188, "And Usage home widget should contain text \"Grok Spawner\"", () => shouldContainText(page, el("Usage home widget"), "Grok Spawner"));
      await session.step(189, "And Usage home widget should contain text \"Grok Connect\"", () => shouldContainText(page, el("Usage home widget"), "Grok Connect"));
      await session.step(190, "And there should be 2 visible line chart viewer in Usage home widget", () => visibleCount(page, 2, el("line chart viewer in Usage home widget")));
      await session.step(191, "And the \"lines\" reading of first line chart viewer in Usage home widget should be 1", () => readingIs(page, "lines", el("first line chart viewer in Usage home widget"), 1));
      await session.step(192, "And the \"rows shown\" reading of first line chart viewer in Usage home widget should be at least 1", () => readingAtLeast(page, "rows shown", el("first line chart viewer in Usage home widget"), 1));
      await session.step(193, "And the \"lines\" reading of second line chart viewer in Usage home widget should be 1", () => readingIs(page, "lines", el("second line chart viewer in Usage home widget"), 1));
      await session.step(194, "And the \"rows shown\" reading of second line chart viewer in Usage home widget should be at least 1", () => readingAtLeast(page, "rows shown", el("second line chart viewer in Usage home widget"), 1));
      await session.step(195, "And \"Open Usage Analysis\" link in Usage home widget should be visible", () => shouldBe(page, el("\"Open Usage Analysis\" link in Usage home widget"), "visible"));
      await session.step(196, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open Usage Analysis in the Usage widget opens the app's Overview view, with nothing logged", async () => {
      await session.step(199, "When user clicks on \"Open Usage Analysis\" link in Usage home widget", () => clickOn(page, el("\"Open Usage Analysis\" link in Usage home widget")));
      await session.step(200, "Then the \"Overview\" view should be open", () => viewOpen(page, "Overview"));
      await session.step(201, "And the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(202, "And no errors should have been logged", () => noErrors(page));
      await session.step(203, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The Overview view Open Usage Analysis opened is brought to the front", async () => {
      await session.step(209, "Then the \"Overview\" view should be current", () => viewIsCurrent(page, "Overview"));
    }, {knownFailure: true});
    await run.scenario("Reports lists recent reports and opens the reports view", async () => {
      await session.step(212, "Then \"Open Reports\" link in Reports home widget should be visible", () => shouldBe(page, el("\"Open Reports\" link in Reports home widget"), "visible"));
      await session.step(213, "When user clicks on \"Open Reports\" link in Reports home widget", () => clickOn(page, el("\"Open Reports\" link in Reports home widget")));
      await session.step(214, "Then the \"Reports\" view should be current", () => viewIsCurrent(page, "Reports"));
      await session.step(215, "And no errors should have been logged", () => noErrors(page));
      await session.step(216, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
