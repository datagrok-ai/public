/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-dock-layout-span.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#239]
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
import {clickOn, enterInto, pressKey, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, clickViewTab, closeAllViews, dialogCloses, dockViewBeside, doubleClickOpensProject, openDataset, openSaveDialog, ownProjectGone, projectsOnServer, saveFromDialog, simpleModeOff, viewDockedBeside, viewSpansPair} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A view docked under two views side by side reopens under those two, not under the side panels", () => {
  const session = feature(test, "features/projects/projects-dock-layout-span.feature", import.meta.url);
  test("A view docked under two views side by side reopens under those two, not under the side panels", {tag: ["@journey", "@realizes:views.projects", "@realizes:core#239"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the user's own project \"BDDLayoutSpan{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDDLayoutSpan{time}")));
    await run.scenario("demog-1000 docked below cars, then iris docked to the right of cars", async () => {
      await session.step(18, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(19, "And user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(20, "And user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(22, "When user drags the tab of the \"demog-1000\" view to the bottom side of the \"cars\" view", () => dockViewBeside(page, "demog-1000", "bottom", "cars"));
      await session.step(23, "Then the \"demog-1000\" view should be docked below the \"cars\" view", () => viewDockedBeside(page, "demog-1000", "below", "cars"));
      await session.step(24, "When user drags the tab of the \"iris\" view to the right side of the \"cars\" view", () => dockViewBeside(page, "iris", "right", "cars"));
      await session.step(25, "Then the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
      await session.step(26, "And the \"demog-1000\" view should span the \"cars\" and \"iris\" views", () => viewSpansPair(page, "demog-1000", "cars", "iris"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reopened dashboard keeps demog-1000 under cars and iris only", async () => {
      await session.step(30, "When user clicks on the tab of the \"cars\" view", () => clickViewTab(page, "cars"));
      await session.step(31, "And user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(32, "And user types \"BDDLayoutSpan{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDDLayoutSpan{time}"), el("text input in \"Save project\" dialog")));
      await session.step(33, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(34, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(35, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(36, "Then 1 project named \"BDDLayoutSpan{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLayoutSpan{time}")));
      await session.step(37, "When user closes all views", () => closeAllViews(page));
      await session.step(38, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(39, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(40, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(41, "And user enters \"BDDLayoutSpan{time}\" into gallery search", () => enterInto(page, session.text("BDDLayoutSpan{time}"), el("gallery search")));
      await session.step(42, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(43, "And user double-clicks on BDDLayoutSpan{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLayoutSpan{time} gallery card"))));
      await session.step(44, "Then the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
      await session.step(45, "And the \"demog-1000\" view should be docked below the \"cars\" view", () => viewDockedBeside(page, "demog-1000", "below", "cars"));
      await session.step(46, "And the \"demog-1000\" view should span the \"cars\" and \"iris\" views", () => viewSpansPair(page, "demog-1000", "cars", "iris"));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
