/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-url-parameters.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20929]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noErrorsButPreviewNoise, onlyRowReads} from '../../bindings/projects-regressions.js';
import {datagrokQuery} from '../../bindings/projects-sources.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, noQueryOnServer, projectsOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: the URL parameters icon of a just-saved dashboard", () => {
  const session = feature(test, "features/projects/projects-regressions-url-parameters.feature", import.meta.url);
  test("Right after the save, Toolbox > Source offers the URL parameters icon", {tag: ["@serial", "@realizes:views.projects", "@known-failure", "@realizes:GROK-20929"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And no query named \"BDDRegIconQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegIconQ{time}")));
    await session.step(24, "And no project named \"BDDRegIcon{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegIcon{time}")));
    await session.step(25, "And a query \"BDDRegIconQ{time}\" on the Datagrok connection is:", () => datagrokQuery(page, session.text("BDDRegIconQ{time}"), "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
    await session.step(30, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(31, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(33, "Given Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(34, "When user double-clicks Databases---Postgres---Datagrok---BDDRegIconQ{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRegIconQ{time} tree node inside browse tree"))));
    await session.step(35, "Then the \"BDDRegIconQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegIconQ{time}")));
    await session.step(36, "And the only row of the table should read \"Project\" in the \"name\" column", () => onlyRowReads(page, "Project", "name"));
    await session.step(37, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(38, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(39, "When user enters \"BDDRegIcon{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegIcon{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(40, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(41, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(42, "And \"Share BDDRegIcon{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegIcon{time}\" dialog")), "visible"));
    await session.step(43, "When user clicks on CANCEL button in \"Share BDDRegIcon{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegIcon{time}\" dialog"))));
    await session.step(44, "Then the \"Share BDDRegIcon{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegIcon{time}")));
    await session.step(45, "And 1 project named \"BDDRegIcon{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDRegIcon{time}")));
    await session.step(46, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    await session.step(47, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(48, "Then REFRESH button inside toolbox should be visible", () => shouldBe(page, el("REFRESH button inside toolbox"), "visible"));
    await knownFailure(async () => {
      await session.step(52, "Then url parameters icon should be visible", () => shouldBe(page, el("url parameters icon"), "visible"));
    });
  });
  test("Reopened from the gallery, Toolbox > Source offers the URL parameters icon", {tag: ["@serial", "@realizes:views.projects", "@realizes:GROK-20929"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And no query named \"BDDRegIconQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegIconQ{time}")));
    await session.step(24, "And no project named \"BDDRegIcon{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegIcon{time}")));
    await session.step(25, "And a query \"BDDRegIconQ{time}\" on the Datagrok connection is:", () => datagrokQuery(page, session.text("BDDRegIconQ{time}"), "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
    await session.step(30, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(31, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(32, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(33, "Given Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(34, "When user double-clicks Databases---Postgres---Datagrok---BDDRegIconQ{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRegIconQ{time} tree node inside browse tree"))));
    await session.step(35, "Then the \"BDDRegIconQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegIconQ{time}")));
    await session.step(36, "And the only row of the table should read \"Project\" in the \"name\" column", () => onlyRowReads(page, "Project", "name"));
    await session.step(37, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(38, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(39, "When user enters \"BDDRegIcon{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegIcon{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(40, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(41, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(42, "And \"Share BDDRegIcon{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegIcon{time}\" dialog")), "visible"));
    await session.step(43, "When user clicks on CANCEL button in \"Share BDDRegIcon{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegIcon{time}\" dialog"))));
    await session.step(44, "Then the \"Share BDDRegIcon{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegIcon{time}")));
    await session.step(45, "And 1 project named \"BDDRegIcon{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDRegIcon{time}")));
    await session.step(46, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    await session.step(47, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(48, "Then REFRESH button inside toolbox should be visible", () => shouldBe(page, el("REFRESH button inside toolbox"), "visible"));
    await session.step(56, "When user closes all views", () => closeAllViews(page));
    await session.step(57, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(58, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(59, "And user enters \"BDDRegIcon{time}\" into gallery search", () => enterInto(page, session.text("BDDRegIcon{time}"), el("gallery search")));
    await session.step(60, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(61, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(62, "When user double-clicks on BDDRegIcon{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegIcon{time} gallery card"))));
    await session.step(63, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
    await session.step(64, "And the \"BDDRegIconQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegIconQ{time}")));
    await session.step(65, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(66, "Then url parameters icon should be visible", () => shouldBe(page, el("url parameters icon"), "visible"));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
});
