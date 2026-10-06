/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/tabs-order-and-projects.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {dragTo, pressKey, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, closeCurrentView, dialogCloses, openDataset, openProjectWithTable, openSaveDialog, ownProjectGone, projectsOnServer, saveFromDialog, simpleModeOff, viewTabsInOrder} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, noTableLeft, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Dataset tabs reordered by drag and drop keep their order through a project", () => {
  const session = feature(test, "features/general/tabs-order-and-projects.feature", import.meta.url);
  test("Dataset tabs reordered by drag and drop keep their order through a project", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(22, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(23, "Then the open table views should be exactly \"browse-import\"", () => openTableViewsExactly(page, "browse-import"));
    await session.step(24, "When user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(25, "Then the open table views should be exactly \"browse-import, cars-small\"", () => openTableViewsExactly(page, "browse-import, cars-small"));
    await session.step(26, "When user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(27, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(28, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await run.scenario("Four datasets, two of them dropped, open as four tabs in that order", async () => {
      await session.step(31, "Then the view tabs should be in the order \"browse-import, cars-small, smiles, curves\"", () => viewTabsInOrder(page, "browse-import, cars-small, smiles, curves"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
      await session.step(33, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A tab dragged onto another moves next to it, and back again", async () => {
      await session.step(36, "When user drags \"curves\" view tab to \"browse-import\" view tab", () => dragTo(page, el("\"curves\" view tab"), el("\"browse-import\" view tab")));
      await session.step(37, "Then the view tabs should be in the order \"browse-import, curves, cars-small, smiles\"", () => viewTabsInOrder(page, "browse-import, curves, cars-small, smiles"));
      await session.step(38, "When user drags \"curves\" view tab to \"smiles\" view tab", () => dragTo(page, el("\"curves\" view tab"), el("\"smiles\" view tab")));
      await session.step(39, "Then the view tabs should be in the order \"browse-import, cars-small, smiles, curves\"", () => viewTabsInOrder(page, "browse-import, cars-small, smiles, curves"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
      await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A dataset opened later goes to the end and leaves the new order alone", async () => {
      await session.step(44, "When user drags \"smiles\" view tab to \"browse-import\" view tab", () => dragTo(page, el("\"smiles\" view tab"), el("\"browse-import\" view tab")));
      await session.step(45, "Then the view tabs should be in the order \"browse-import, smiles, cars-small, curves\"", () => viewTabsInOrder(page, "browse-import, smiles, cars-small, curves"));
      await session.step(46, "When user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(47, "Then the view tabs should be in the order \"browse-import, smiles, cars-small, curves, demog\"", () => viewTabsInOrder(page, "browse-import, smiles, cars-small, curves, demog"));
      await session.step(48, "When user closes the current view", () => closeCurrentView(page));
      await session.step(49, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
      await session.step(50, "And the view tabs should be in the order \"browse-import, smiles, cars-small, curves\"", () => viewTabsInOrder(page, "browse-import, smiles, cars-small, curves"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The reordered tabs come back in that order from the saved project", async () => {
      await session.step(54, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(55, "And user types \"BDD-Tabs-{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDD-Tabs-{time}"), el("text input in \"Save project\" dialog")));
      await session.step(56, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(57, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(58, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(59, "Then 1 project named \"BDD-Tabs-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Tabs-{time}")));
      await session.step(60, "When user closes all views", () => closeAllViews(page));
      await session.step(61, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(62, "When user opens the \"BDD-Tabs-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Tabs-{time}")));
      await session.step(63, "Then the view tabs should be in the order \"browse-import, smiles, cars-small, curves\"", () => viewTabsInOrder(page, "browse-import, smiles, cars-small, curves"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
