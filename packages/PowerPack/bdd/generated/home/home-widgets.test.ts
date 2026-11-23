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
import '../../bindings/identifiers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {adminAccountOnServer, everyWidgetHasContent, homeWidgetsAre, tipOfTheDay, widgetStoredAs, widgetsAllShown} from '../../bindings/home.js';
import {searchFinished} from '../../bindings/search.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, containsOneOf, followingShouldBe, hoverOver, placeholderStarts, shouldBe, shouldContainText, shouldHaveText, shouldNotBe, tabShowing, typeInto, uncheck, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, openDataset, signInAs, switchView, urlShouldContain, urlShouldNotContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, readingAtLeast, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The widgets of the Home page", () => {
  const session = feature(test, "features/home/home-widgets.feature", import.meta.url);
  test("The widgets of the Home page", {tag: ["@journey", "@serial", "@realizes:powerpack.view.welcome", "@realizes:powerpack.dashboard.spotlight", "@realizes:powerpack.dashboard.community"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(50, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "And an administrator account \"bddhomeadmin\" is on the server", () => adminAccountOnServer(page, "bddhomeadmin"));
    await session.step(52, "And user signs in as \"bddhomeadmin\"", () => signInAs(page, "bddhomeadmin"));
    await session.step(53, "And every widget of the Home page is stored as shown, now and when the feature ends", () => widgetsAllShown(page));
    await run.scenario("The Home page shows the search box and four loaded widgets", async () => {
      await session.step(56, "Then home search should be visible", () => shouldBe(page, el("home search"), "visible"));
      await session.step(57, "And the placeholder of home search should start with \"Search everywhere\"", () => placeholderStarts(page, el("home search"), "Search everywhere"));
      await session.step(58, "And the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(59, "And every widget of the Home page should show content", () => everyWidgetHasContent(page));
      await session.step(60, "And title of Community home widget should have text \"Community\"", () => shouldHaveText(page, el("title of Community home widget"), "Community"));
      await session.step(61, "And title of Spotlight home widget should be hidden", () => shouldBe(page, el("title of Spotlight home widget"), "hidden"));
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The close icon shows on hover and removes the widget", async () => {
      await session.step(66, "Then close icon of Community home widget should be present", () => shouldBe(page, el("close icon of Community home widget"), "present"));
      await session.step(67, "And close icon of Community home widget should be hidden", () => shouldBe(page, el("close icon of Community home widget"), "hidden"));
      await session.step(68, "When user hovers over Community home widget", () => hoverOver(page, el("Community home widget")));
      await session.step(69, "Then close icon of Community home widget should be visible", () => shouldBe(page, el("close icon of Community home widget"), "visible"));
      await session.step(70, "When user hovers over close icon of Community home widget", () => hoverOver(page, el("close icon of Community home widget")));
      await session.step(71, "Then tooltip should contain text \"Remove\"", () => shouldContainText(page, el("tooltip"), "Remove"));
      await session.step(72, "When user clicks on close icon of Community home widget", () => clickOn(page, el("close icon of Community home widget")));
      await session.step(73, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(74, "And the Home page should show the widgets \"Spotlight, Reports, Usage\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage"));
      await session.step(75, "And the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(76, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(77, "And user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(78, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(79, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Customize form hides and shows a widget", async () => {
      await session.step(83, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(84, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Spotlight\" input in context panel"],["\"Community\" input in context panel"],["\"Usage\" input in context panel"],["\"Reports\" input in context panel"]]), [["\"Spotlight\" input in context panel"],["\"Community\" input in context panel"],["\"Usage\" input in context panel"],["\"Reports\" input in context panel"]]);
      await session.step(89, "And \"Community\" input in context panel should be checked", () => shouldBe(page, el("\"Community\" input in context panel"), "checked"));
      await session.step(90, "When user unchecks \"Community\" input in context panel", () => uncheck(page, el("\"Community\" input in context panel")));
      await session.step(91, "Then Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(92, "And the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(93, "When user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(94, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(95, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A hidden widget stays hidden after a reload", async () => {
      await session.step(99, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(100, "And user unchecks \"Community\" input in context panel", () => uncheck(page, el("\"Community\" input in context panel")));
      await session.step(101, "Then the Community widget should be stored as hidden", () => widgetStoredAs(page, "Community", "hidden"));
      await session.step(102, "When user reloads the page", () => reloadPage(page));
      await session.step(103, "Then the Home page should show the widgets \"Spotlight, Reports, Usage\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage"));
      await session.step(104, "And Community home widget should be absent", () => shouldBe(page, el("Community home widget"), "absent"));
      await session.step(105, "When user clicks on \"Customize widgets...\" link", () => clickOn(page, el("\"Customize widgets...\" link")));
      await session.step(106, "Then \"Community\" input in context panel should not be checked", () => shouldNotBe(page, el("\"Community\" input in context panel"), "checked"));
      await session.step(107, "When user checks \"Community\" input in context panel", () => check(page, el("\"Community\" input in context panel")));
      await session.step(108, "Then Community home widget should be visible", () => shouldBe(page, el("Community home widget"), "visible"));
      await session.step(109, "And the Community widget should be stored as shown", () => widgetStoredAs(page, "Community", "shown"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A search replaces the widgets, and clearing it brings them back", async () => {
      await session.step(113, "Then home widgets panel should be visible", () => shouldBe(page, el("home widgets panel"), "visible"));
      await session.step(114, "When user types \"demog\" into home search", () => typeInto(page, "demog", el("home search")));
      await session.step(115, "Then home widgets panel should be hidden", () => shouldBe(page, el("home widgets panel"), "hidden"));
      await session.step(116, "And home search results should be visible", () => shouldBe(page, el("home search results"), "visible"));
      await session.step(117, "And the page address should contain \"search?q=demog\"", () => urlShouldContain(page, "search?q=demog"));
      await session.step(118, "And the search should have finished", () => searchFinished(page));
      await session.step(119, "When user clears home search", () => clearField(page, el("home search")));
      await session.step(120, "Then home widgets panel should be visible", () => shouldBe(page, el("home widgets panel"), "visible"));
      await session.step(121, "And home search results should be hidden", () => shouldBe(page, el("home search results"), "hidden"));
      await session.step(122, "And the page address should not contain \"?q=\"", () => urlShouldNotContain(page, "?q="));
      await session.step(123, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Home icon brings back the same widgets", async () => {
      await session.step(126, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(127, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(128, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(129, "When user clicks on \"Home\" icon in browse toolbar", () => clickOn(page, el("\"Home\" icon in browse toolbar")));
      await session.step(130, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(131, "And the Home page should show the widgets \"Spotlight, Reports, Usage, Community\"", () => homeWidgetsAre(page, "Spotlight, Reports, Usage, Community"));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Spotlight has six tabs, and each shows its own page", async () => {
      await session.step(135, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Workspace tab in Spotlight home widget"],["Spotlight tab in Spotlight home widget"],["Favorites tab in Spotlight home widget"],["Notifications tab in Spotlight home widget"],["\"My Activity\" tab in Spotlight home widget"],["Learn tab in Spotlight home widget"]]), [["Workspace tab in Spotlight home widget"],["Spotlight tab in Spotlight home widget"],["Favorites tab in Spotlight home widget"],["Notifications tab in Spotlight home widget"],["\"My Activity\" tab in Spotlight home widget"],["Learn tab in Spotlight home widget"]]);
      await session.step(142, "When user clicks on Spotlight tab in Spotlight home widget", () => clickOn(page, el("Spotlight tab in Spotlight home widget")));
      await session.step(143, "Then the \"Spotlight\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Spotlight", el("Spotlight home widget")));
      await session.step(144, "And spotlight page of Spotlight home widget should contain one of the texts \"Action required | Shared with me | Recent | Interactive Tutorials\"", () => containsOneOf(page, el("spotlight page of Spotlight home widget"), "Action required | Shared with me | Recent | Interactive Tutorials"));
      await session.step(145, "When user clicks on Favorites tab in Spotlight home widget", () => clickOn(page, el("Favorites tab in Spotlight home widget")));
      await session.step(146, "Then the \"Favorites\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Favorites", el("Spotlight home widget")));
      await session.step(147, "When user clicks on Notifications tab in Spotlight home widget", () => clickOn(page, el("Notifications tab in Spotlight home widget")));
      await session.step(148, "Then the \"Notifications\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Notifications", el("Spotlight home widget")));
      await session.step(149, "When user clicks on \"My Activity\" tab in Spotlight home widget", () => clickOn(page, el("\"My Activity\" tab in Spotlight home widget")));
      await session.step(150, "Then the \"My Activity\" tab of Spotlight home widget should be showing", () => tabShowing(page, "My Activity", el("Spotlight home widget")));
      await session.step(151, "When user clicks on Learn tab in Spotlight home widget", () => clickOn(page, el("Learn tab in Spotlight home widget")));
      await session.step(152, "Then the \"Learn\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Learn", el("Spotlight home widget")));
      await session.step(153, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["VIDEO tab in Spotlight home widget"],["WIKI tab in Spotlight home widget"],["DEMO tab in Spotlight home widget"],["TUTORIALS tab in Spotlight home widget"]]), [["VIDEO tab in Spotlight home widget"],["WIKI tab in Spotlight home widget"],["DEMO tab in Spotlight home widget"],["TUTORIALS tab in Spotlight home widget"]]);
      await session.step(158, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Meetings\" text in Spotlight home widget"],["\"Develop\" text in Spotlight home widget"],["\"Cheminformatics\" text in Spotlight home widget"],["\"Visualize\" text in Spotlight home widget"],["\"Explore\" text in Spotlight home widget"]]), [["\"Meetings\" text in Spotlight home widget"],["\"Develop\" text in Spotlight home widget"],["\"Cheminformatics\" text in Spotlight home widget"],["\"Visualize\" text in Spotlight home widget"],["\"Explore\" text in Spotlight home widget"]]);
      await session.step(164, "When user clicks on WIKI tab in Spotlight home widget", () => clickOn(page, el("WIKI tab in Spotlight home widget")));
      await session.step(165, "Then the \"WIKI\" tab of Spotlight home widget should be showing", () => tabShowing(page, "WIKI", el("Spotlight home widget")));
      await session.step(166, "And the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Bioinformatics\" text in Spotlight home widget"],["\"Overview\" text in Spotlight home widget"],["\"Access\" text in Spotlight home widget"],["\"Transform\" text in Spotlight home widget"],["\"Compute\" text in Spotlight home widget"],["\"Govern\" text in Spotlight home widget"]]), [["\"Bioinformatics\" text in Spotlight home widget"],["\"Overview\" text in Spotlight home widget"],["\"Access\" text in Spotlight home widget"],["\"Transform\" text in Spotlight home widget"],["\"Compute\" text in Spotlight home widget"],["\"Govern\" text in Spotlight home widget"]]);
      await session.step(173, "When user clicks on DEMO tab in Spotlight home widget", () => clickOn(page, el("DEMO tab in Spotlight home widget")));
      await session.step(174, "Then the \"DEMO\" tab of Spotlight home widget should be showing", () => tabShowing(page, "DEMO", el("Spotlight home widget")));
      await session.step(175, "When user clicks on TUTORIALS tab in Spotlight home widget", () => clickOn(page, el("TUTORIALS tab in Spotlight home widget")));
      await session.step(176, "Then the \"TUTORIALS\" tab of Spotlight home widget should be showing", () => tabShowing(page, "TUTORIALS", el("Spotlight home widget")));
      await session.step(177, "When user clicks on VIDEO tab in Spotlight home widget", () => clickOn(page, el("VIDEO tab in Spotlight home widget")));
      await session.step(178, "Then the \"VIDEO\" tab of Spotlight home widget should be showing", () => tabShowing(page, "VIDEO", el("Spotlight home widget")));
      await session.step(179, "When user clicks on Workspace tab in Spotlight home widget", () => clickOn(page, el("Workspace tab in Spotlight home widget")));
      await session.step(180, "Then the \"Workspace\" tab of Spotlight home widget should be showing", () => tabShowing(page, "Workspace", el("Spotlight home widget")));
      await session.step(181, "And no errors should have been logged", () => noErrors(page));
      await session.step(182, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The tip of the day at the bottom of Spotlight opens what it names", async () => {
      await session.step(185, "Then the tip of the day of Spotlight home widget should open what it names", () => tipOfTheDay(page, el("Spotlight home widget")));
      await session.step(186, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Usage shows its user and error charts and the state of the services", async () => {
      await session.step(189, "Then Usage home widget should contain text \"Users\"", () => shouldContainText(page, el("Usage home widget"), "Users"));
      await session.step(190, "And Usage home widget should contain text \"Errors\"", () => shouldContainText(page, el("Usage home widget"), "Errors"));
      await session.step(191, "And Usage home widget should contain text \"System\"", () => shouldContainText(page, el("Usage home widget"), "System"));
      await session.step(192, "And there should be 2 visible line chart viewer in Usage home widget", () => visibleCount(page, 2, el("line chart viewer in Usage home widget")));
      await session.step(193, "And the \"lines\" reading of first line chart viewer in Usage home widget should be 1", () => readingIs(page, "lines", el("first line chart viewer in Usage home widget"), 1));
      await session.step(194, "And the \"rows shown\" reading of first line chart viewer in Usage home widget should be at least 1", () => readingAtLeast(page, "rows shown", el("first line chart viewer in Usage home widget"), 1));
      await session.step(195, "And the \"lines\" reading of second line chart viewer in Usage home widget should be 1", () => readingIs(page, "lines", el("second line chart viewer in Usage home widget"), 1));
      await session.step(196, "And the \"rows shown\" reading of second line chart viewer in Usage home widget should be at least 1", () => readingAtLeast(page, "rows shown", el("second line chart viewer in Usage home widget"), 1));
      await session.step(197, "And \"Open Usage Analysis\" link in Usage home widget should be visible", () => shouldBe(page, el("\"Open Usage Analysis\" link in Usage home widget"), "visible"));
      await session.step(198, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open Usage Analysis in the Usage widget opens the app's Overview view in front, with nothing logged", async () => {
      await session.step(201, "When user clicks on \"Open Usage Analysis\" link in Usage home widget", () => clickOn(page, el("\"Open Usage Analysis\" link in Usage home widget")));
      await session.step(202, "Then the \"Overview\" view should be current", () => viewIsCurrent(page, "Overview"));
      await session.step(203, "And no errors should have been logged", () => noErrors(page));
      await session.step(204, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(205, "When user switches to the \"Home\" view", () => switchView(page, "Home"));
      await session.step(206, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    });
    await run.scenario("The Reports widget opens the reports view", async () => {
      await session.step(209, "Then \"Open Reports\" link in Reports home widget should be visible", () => shouldBe(page, el("\"Open Reports\" link in Reports home widget"), "visible"));
      await session.step(210, "When user clicks on \"Open Reports\" link in Reports home widget", () => clickOn(page, el("\"Open Reports\" link in Reports home widget")));
      await session.step(211, "Then the \"Reports\" view should be current", () => viewIsCurrent(page, "Reports"));
      await session.step(212, "And no errors should have been logged", () => noErrors(page));
      await session.step(213, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
