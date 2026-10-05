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
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {dragTo, pressKey, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, openDataset, openProjectWithTable, openSaveDialog, ownProjectGone, projectsOnServer, saveFromDialog, savedAsSnapshot, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, noTableLeft, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Dataset tabs reordered by drag and drop keep their order through a project", () => {
  const session = feature(test, "features/general/tabs-order-and-projects.feature", import.meta.url);
  test("Four datasets, two of them dropped, open as four tabs in that order", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(24, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(25, "Then the open table views should be exactly \"browse-import\"", () => openTableViewsExactly(page, "browse-import"));
    await session.step(26, "When user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(27, "Then the open table views should be exactly \"browse-import, cars-small\"", () => openTableViewsExactly(page, "browse-import, cars-small"));
    await session.step(28, "When user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(29, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(30, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(33, "Then 4th view tab should have text \"browse-import\"", () => shouldHaveText(page, el("4th view tab"), "browse-import"));
    await session.step(34, "And 5th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("5th view tab"), "cars-small"));
    await session.step(35, "And 6th view tab should have text \"smiles\"", () => shouldHaveText(page, el("6th view tab"), "smiles"));
    await session.step(36, "And 7th view tab should have text \"curves\"", () => shouldHaveText(page, el("7th view tab"), "curves"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
    await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A tab dragged onto another moves next to it, and back again", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(24, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(25, "Then the open table views should be exactly \"browse-import\"", () => openTableViewsExactly(page, "browse-import"));
    await session.step(26, "When user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(27, "Then the open table views should be exactly \"browse-import, cars-small\"", () => openTableViewsExactly(page, "browse-import, cars-small"));
    await session.step(28, "When user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(29, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(30, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(41, "When user drags \"curves\" view tab to \"browse-import\" view tab", () => dragTo(page, el("\"curves\" view tab"), el("\"browse-import\" view tab")));
    await session.step(42, "Then 4th view tab should have text \"browse-import\"", () => shouldHaveText(page, el("4th view tab"), "browse-import"));
    await session.step(43, "And 5th view tab should have text \"curves\"", () => shouldHaveText(page, el("5th view tab"), "curves"));
    await session.step(44, "And 6th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("6th view tab"), "cars-small"));
    await session.step(45, "And 7th view tab should have text \"smiles\"", () => shouldHaveText(page, el("7th view tab"), "smiles"));
    await session.step(46, "When user drags \"curves\" view tab to \"smiles\" view tab", () => dragTo(page, el("\"curves\" view tab"), el("\"smiles\" view tab")));
    await session.step(47, "Then 4th view tab should have text \"browse-import\"", () => shouldHaveText(page, el("4th view tab"), "browse-import"));
    await session.step(48, "And 5th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("5th view tab"), "cars-small"));
    await session.step(49, "And 6th view tab should have text \"smiles\"", () => shouldHaveText(page, el("6th view tab"), "smiles"));
    await session.step(50, "And 7th view tab should have text \"curves\"", () => shouldHaveText(page, el("7th view tab"), "curves"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
    await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A dataset opened later goes to the end and leaves the new order alone", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(24, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(25, "Then the open table views should be exactly \"browse-import\"", () => openTableViewsExactly(page, "browse-import"));
    await session.step(26, "When user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(27, "Then the open table views should be exactly \"browse-import, cars-small\"", () => openTableViewsExactly(page, "browse-import, cars-small"));
    await session.step(28, "When user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(29, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(30, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(55, "When user drags \"smiles\" view tab to \"browse-import\" view tab", () => dragTo(page, el("\"smiles\" view tab"), el("\"browse-import\" view tab")));
    await session.step(56, "Then 5th view tab should have text \"smiles\"", () => shouldHaveText(page, el("5th view tab"), "smiles"));
    await session.step(57, "And 6th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("6th view tab"), "cars-small"));
    await session.step(58, "When user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(59, "Then 4th view tab should have text \"browse-import\"", () => shouldHaveText(page, el("4th view tab"), "browse-import"));
    await session.step(60, "And 5th view tab should have text \"smiles\"", () => shouldHaveText(page, el("5th view tab"), "smiles"));
    await session.step(61, "And 6th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("6th view tab"), "cars-small"));
    await session.step(62, "And 7th view tab should have text \"curves\"", () => shouldHaveText(page, el("7th view tab"), "curves"));
    await session.step(63, "And 8th view tab should have text \"demog\"", () => shouldHaveText(page, el("8th view tab"), "demog"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("The reordered tabs come back in that order from the saved project", {tag: ["@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(23, "And the user's own project \"BDD-Tabs-{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDD-Tabs-{time}")));
    await session.step(24, "When user drops the \"fixtures/browse-import.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/browse-import.csv", el("status bar")));
    await session.step(25, "Then the open table views should be exactly \"browse-import\"", () => openTableViewsExactly(page, "browse-import"));
    await session.step(26, "When user drops the \"fixtures/cars-small.csv\" file of the project onto status bar", () => dropFile(page, "fixtures/cars-small.csv", el("status bar")));
    await session.step(27, "Then the open table views should be exactly \"browse-import, cars-small\"", () => openTableViewsExactly(page, "browse-import, cars-small"));
    await session.step(28, "When user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await session.step(29, "And user opens curves dataset", () => openDataset(page, ds("curves")));
    await session.step(30, "Then the open table views should be exactly \"browse-import, cars-small, smiles, curves\"", () => openTableViewsExactly(page, "browse-import, cars-small, smiles, curves"));
    await session.step(67, "When user drags \"smiles\" view tab to \"browse-import\" view tab", () => dragTo(page, el("\"smiles\" view tab"), el("\"browse-import\" view tab")));
    await session.step(68, "Then 5th view tab should have text \"smiles\"", () => shouldHaveText(page, el("5th view tab"), "smiles"));
    await session.step(69, "And 6th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("6th view tab"), "cars-small"));
    await session.step(70, "When user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
    await session.step(71, "And user types \"BDD-Tabs-{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDD-Tabs-{time}"), el("text input in \"Save project\" dialog")));
    await session.step(72, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(73, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(74, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(75, "Then 1 project named \"BDD-Tabs-{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDD-Tabs-{time}")));
    await session.step(76, "And the \"browse-import\" table of the \"BDD-Tabs-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "browse-import", session.text("BDD-Tabs-{time}")));
    await session.step(77, "And the \"cars-small\" table of the \"BDD-Tabs-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "cars-small", session.text("BDD-Tabs-{time}")));
    await session.step(78, "And the \"smiles\" table of the \"BDD-Tabs-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "smiles", session.text("BDD-Tabs-{time}")));
    await session.step(79, "And the \"curves\" table of the \"BDD-Tabs-{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "curves", session.text("BDD-Tabs-{time}")));
    await session.step(80, "When user closes all views", () => closeAllViews(page));
    await session.step(81, "Then no table should be left in the workspace", () => noTableLeft(page));
    await session.step(82, "When user reloads the page", () => reloadPage(page));
    await session.step(83, "And user opens the \"BDD-Tabs-{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDD-Tabs-{time}")));
    await session.step(84, "Then 4th view tab should have text \"browse-import\"", () => shouldHaveText(page, el("4th view tab"), "browse-import"));
    await session.step(85, "And 5th view tab should have text \"smiles\"", () => shouldHaveText(page, el("5th view tab"), "smiles"));
    await session.step(86, "And 6th view tab should have text \"cars-small\"", () => shouldHaveText(page, el("6th view tab"), "cars-small"));
    await session.step(87, "And 7th view tab should have text \"curves\"", () => shouldHaveText(page, el("7th view tab"), "curves"));
    await session.step(88, "And no errors should have been logged", () => noErrors(page));
  });
});
