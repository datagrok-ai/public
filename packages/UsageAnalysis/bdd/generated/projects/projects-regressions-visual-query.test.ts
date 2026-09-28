/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-visual-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20535]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clearSavedParameters} from '../../bindings/pivot-table.js';
import {noErrorsButPreviewNoise, onlyRowReads} from '../../bindings/projects-regressions.js';
import {addToBuilderRow, builderRowCount, builderRowHolds, cleanQueryLayout, savedQueryFilters} from '../../bindings/queries.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, noProjectOnServer, noQueryOnServer, queriesOnServer, reloadedByDataSync, savedWithDataSync, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a project of a visual query with a parameter", () => {
  const session = feature(test, "features/projects/projects-regressions-visual-query.feature", import.meta.url);
  test("The reopened visual query project holds the row of the value it was run with", {tag: ["@serial", "@realizes:views.projects", "@known-failure", "@realizes:GROK-20535"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(32, "And no query named \"BDDRegVQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegVQ{time}")));
    await session.step(33, "And no project named \"BDDRegVQProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegVQProj{time}")));
    await session.step(34, "And the layout saved for the query \"BDDRegVQ{time}\" is deleted at the end", () => cleanQueryLayout(page, session.text("BDDRegVQ{time}")));
    await session.step(35, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(36, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(37, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(38, "And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas tree node inside browse tree")));
    await session.step(39, "And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree")));
    await session.step(40, "When user picks \"New Visual Query...\" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity-types tree node inside browse tree", () => pickFromContextMenu(page, "New Visual Query...", el("Databases---Postgres---Datagrok---Schemas---public---entity-types tree node inside browse tree")));
    await session.step(41, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(42, "When user adds \"name\" to the \"Where\" row of the visual query", () => addToBuilderRow(page, "name", "Where"));
    await session.step(43, "And user enters \"Project\" into visual query name condition", () => enterInto(page, "Project", el("visual query name condition")));
    await session.step(44, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(45, "And user checks visual query name parameter checkbox", () => check(page, el("visual query name parameter checkbox")));
    await session.step(46, "Then visual query name parameter checkbox should be checked", () => shouldBe(page, el("visual query name parameter checkbox"), "checked"));
    await session.step(47, "When user adds \"name\" to the \"Order by\" row of the visual query", () => addToBuilderRow(page, "name", "Order by"));
    await session.step(48, "Then the \"Order by\" row of the visual query should hold \"name\"", () => builderRowHolds(page, "Order by", "name"));
    await session.step(49, "And the visual query should have run to 1 row", () => builderRowCount(page, 1));
    await session.step(50, "When user enters \"BDDRegVQ{time}\" into Name input", () => enterInto(page, session.text("BDDRegVQ{time}"), el("Name input")));
    await session.step(51, "And user clicks on Save button", () => clickOn(page, el("Save button")));
    await session.step(52, "Then 1 query named \"BDDRegVQ{time}\" should be on the server", () => queriesOnServer(page, 1, session.text("BDDRegVQ{time}")));
    await session.step(53, "And the query \"BDDRegVQ{time}\" on the server should filter \"name\" by \"Project\" as a parameter", () => savedQueryFilters(page, session.text("BDDRegVQ{time}"), "name", "Project"));
    await session.step(54, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(55, "When user clicks on \"Run query...\" action in toolbox", () => clickOn(page, el("\"Run query...\" action in toolbox")));
    await session.step(56, "Then \"BDDRegVQ{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDRegVQ{time}\" dialog")), "visible"));
    await session.step(57, "And Name text input in \"BDDRegVQ{time}\" dialog should have value \"Project\"", () => shouldHaveValue(page, el(session.text("Name text input in \"BDDRegVQ{time}\" dialog")), "Project"));
    await session.step(58, "When user enters \"Script\" into Name text input in \"BDDRegVQ{time}\" dialog", () => enterInto(page, "Script", el(session.text("Name text input in \"BDDRegVQ{time}\" dialog"))));
    await session.step(59, "And user clicks on OK button in \"BDDRegVQ{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDDRegVQ{time}\" dialog"))));
    await session.step(60, "Then the \"BDDRegVQ{time}\" dialog should close", () => dialogCloses(page, session.text("BDDRegVQ{time}")));
    await session.step(61, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(62, "And the only row of the table should read \"Script\" in the \"name\" column", () => onlyRowReads(page, "Script", "name"));
    await session.step(63, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(64, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(65, "When user enters \"BDDRegVQProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegVQProj{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(66, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(67, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(68, "And \"Share BDDRegVQProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegVQProj{time}\" dialog")), "visible"));
    await session.step(69, "When user clicks on CANCEL button in \"Share BDDRegVQProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegVQProj{time}\" dialog"))));
    await session.step(70, "Then the \"Share BDDRegVQProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegVQProj{time}")));
    await session.step(71, "And the \"BDDRegVQ{time}\" table of the \"BDDRegVQProj{time}\" project should be saved with data sync", () => savedWithDataSync(page, session.text("BDDRegVQ{time}"), session.text("BDDRegVQProj{time}")));
    await session.step(72, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    await session.step(73, "When user closes all views", () => closeAllViews(page));
    await session.step(74, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(75, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(76, "And user enters \"BDDRegVQProj{time}\" into gallery search", () => enterInto(page, session.text("BDDRegVQProj{time}"), el("gallery search")));
    await session.step(77, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(78, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(79, "When user double-clicks on BDDRegVQProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegVQProj{time} gallery card"))));
    await session.step(80, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
    await session.step(81, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(82, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
    await knownFailure(async () => {
      await session.step(86, "Then the only row of the table should read \"Script\" in the \"name\" column", () => onlyRowReads(page, "Script", "name"));
    });
  });
});
