/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-result-project.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, currentViewType, dialogCloses, noProjectOnServer, openProjectWithTable, projectsOnServer, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A query result saved as a project", () => {
  const session = feature(test, "features/queries/query-result-project.feature", import.meta.url);
  test("A query result saved as a project", {tag: ["@journey", "@serial", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And no project named \"BDD-Q-proj-{time}\" is on the server", () => noProjectOnServer(page, session.text("BDD-Q-proj-{time}")));
    await run.scenario("A query result with a viewer is saved as a project", async () => {
      await session.step(23, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(24, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(25, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
      await session.step(27, "When user hovers over Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => hoverOver(page, el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
      await session.step(28, "And user picks \"Run\" from the context menu of Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---Postgres---NorthwindTest---PostgresByStringChoices tree node inside browse tree")));
      await session.step(29, "And user selects \"France\" in \"Ship Country\" input in \"PostgresByStringChoices\" dialog", () => selectIn(page, "France", el("\"Ship Country\" input in \"PostgresByStringChoices\" dialog")));
      await session.step(30, "And user clicks on OK button in \"PostgresByStringChoices\" dialog", () => clickOn(page, el("OK button in \"PostgresByStringChoices\" dialog")));
      await session.step(31, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(32, "And the table should have 77 rows", () => rowCount(page, 77));
      await session.step(33, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(34, "When user clicks on trellis plot icon in toolbox", () => clickOn(page, el("trellis plot icon in toolbox")));
      await session.step(35, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(36, "When user clicks on Save button", () => clickOn(page, el("Save button")));
      await session.step(37, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(38, "When user enters \"BDD-Q-proj-{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDD-Q-proj-{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(39, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(40, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(41, "And 1 project named \"BDD-Q-proj-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Q-proj-{time}")));
      await session.step(42, "And \"Share BDD-Q-proj-{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDD-Q-proj-{time}\" dialog")), "visible"));
      await session.step(43, "When user clicks on CANCEL button in \"Share BDD-Q-proj-{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDD-Q-proj-{time}\" dialog"))));
      await session.step(44, "Then the \"Share BDD-Q-proj-{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDD-Q-proj-{time}")));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The saved project opens with its table and its viewer", async () => {
      await session.step(49, "When user closes all views", () => closeAllViews(page));
      await session.step(50, "And user opens the \"BDD-Q-proj-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Q-proj-{time}")));
      await session.step(51, "Then the table should have 77 rows", () => rowCount(page, 77));
      await session.step(52, "And the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
