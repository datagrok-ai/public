/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/tabs-order-and-projects.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, dialogCloses, openDataset, openProjectWithTable, openSaveDialog, ownProjectGone, projectsOnServer, saveFromDialog} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, noTableLeft, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Datasets opened together keep their order through a project", () => {
  const session = feature(test, "features/general/tabs-order-and-projects.feature", import.meta.url);
  test("Four datasets, two of them dropped, open in four views without an error", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(22, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(23, "And user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(24, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(25, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(26, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
    await session.step(28, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The project saved through the ribbon opens its views in the same order", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(31, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(32, "And user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(33, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(34, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(35, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(36, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
    await session.step(37, "And user types \"BDD-Tabs-{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDD-Tabs-{time}"), el("text input in \"Save project\" dialog")));
    await session.step(38, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(39, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(40, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(41, "Then 1 project named \"BDD-Tabs-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Tabs-{time}")));
    await session.step(42, "When user closes all views", () => closeAllViews(page));
    await session.step(43, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(44, "When user opens the \"BDD-Tabs-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Tabs-{time}")));
    await session.step(45, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(46, "When user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(47, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves, demog\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves, demog"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
