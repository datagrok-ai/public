/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-quick-save.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#945]
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
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {dblclickRunsQuery} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelOpen, dialogCloses, doubleClickOpensProject, noProjectOnServer, queryTextOnServer, reloadedByDataSync, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: OK clicked as soon as the Save dialog opens", () => {
  const session = feature(test, "features/projects/projects-regressions-quick-save.feature", import.meta.url);
  test("The dashboard saved with an immediate OK reopens with its table", {tag: ["@realizes:views.projects", "@realizes:core#945"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And no project named \"BDDQuickSave{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDQuickSave{time}")));
    await session.step(17, "And a query \"BDDQuickSaveQ{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQuickSaveQ{time}"), "System:Datagrok", "--input: string typeName = \"Script\"\nselect name from entity_types where name = @typeName"));
    await session.step(25, "When user reloads the page", () => reloadPage(page));
    await session.step(26, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(28, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(30, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(31, "When user double-clicks on Databases---Postgres---Datagrok---BDDQuickSaveQ{time} tree node inside browse tree and the query it runs completes", () => dblclickRunsQuery(page, el(session.text("Databases---Postgres---Datagrok---BDDQuickSaveQ{time} tree node inside browse tree"))));
    await session.step(32, "Then the \"BDDQuickSaveQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDQuickSaveQ{time}")));
    await session.step(33, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(34, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(35, "When user enters \"BDDQuickSave{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDQuickSave{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(37, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(38, "And \"Share BDDQuickSave{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDQuickSave{time}\" dialog")), "visible"));
    await session.step(39, "When user clicks on CANCEL button in \"Share BDDQuickSave{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDQuickSave{time}\" dialog"))));
    await session.step(40, "Then the \"Share BDDQuickSave{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDQuickSave{time}")));
    await session.step(41, "When user closes all views", () => closeAllViews(page));
    await session.step(42, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(44, "And user enters \"BDDQuickSave{time}\" into gallery search", () => enterInto(page, session.text("BDDQuickSave{time}"), el("gallery search")));
    await session.step(45, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(46, "And user double-clicks on BDDQuickSave{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDQuickSave{time} gallery card"))));
    await session.step(47, "Then the \"BDDQuickSaveQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDQuickSaveQ{time}")));
    await session.step(48, "And the table should have 1 rows", () => rowCount(page, 1));
    await session.step(49, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
  });
});
