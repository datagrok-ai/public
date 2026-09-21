/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {axesShouldBe} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, legendLists, loadLayout, noErrors, pickFromContextMenu, propertiesShouldBe, propertyShouldBe, propertyShouldNotBe, saveLayoutToServer, setProperties, setProperty, showsRows, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot persistence and Pick Up / Apply", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-persistence.feature", import.meta.url);
  test("PC plot persistence and Pick Up / Apply", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"],["Color","RACE"],["Show Title","true"],["Title","PC Persistence Probe"]]));
    await session.step(23, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(24, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(25, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
    await run.scenario("A saved layout restores the viewer set and the configuration", async () => {
      await session.step(28, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(29, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(30, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(31, "When user loads the saved layout", () => loadLayout(page));
      await session.step(32, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
      await session.step(33, "And the open tableview should have 1 pc plot viewer", () => viewerCount(page, 1, "pc plot"));
      await session.step(34, "And properties of pc plot viewer should be:", () => propertiesShouldBe(page, el("pc plot viewer"), [["Column Names","AGE, HEIGHT, WEIGHT"],["Color","RACE"],["Title","PC Persistence Probe"]]));
      await session.step(38, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(39, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(40, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(41, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A saved project survives Close All and a reopen", async () => {
      await session.step(45, "When user saves the current view as project \"zz-pcplot-bdd-probe\"", () => saveAsProject(page, "zz-pcplot-bdd-probe"));
      await session.step(46, "And user closes all views", () => closeAllViews(page));
      await session.step(47, "And user opens the \"zz-pcplot-bdd-probe\" project", () => openProject(page, "zz-pcplot-bdd-probe"));
      await session.step(48, "Then pc plot viewer should be visible", () => shouldBe(page, el("pc plot viewer"), "visible"));
      await session.step(49, "And properties of pc plot viewer should be:", () => propertiesShouldBe(page, el("pc plot viewer"), [["Column Names","AGE, HEIGHT, WEIGHT"],["Color","RACE"],["Title","PC Persistence Probe"]]));
      await session.step(53, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(54, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(55, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pick Up on one plot and Apply on another copies its settings, once", async () => {
      await session.step(59, "When user adds a pc plot viewer", () => addViewer(page, "pc plot"));
      await session.step(60, "Then the open tableview should have 2 pc plot viewers", () => viewerCount(page, 2, "pc plot"));
      await session.step(61, "When user sets properties of first pc plot viewer:", () => setProperties(page, el("first pc plot viewer"), [["Column Names","AGE, WEIGHT, STARTED"],["Log Columns","AGE"],["Color","RACE"],["Legend Position","Left"],["Title","Source Plot"]]));
      await session.step(67, "Then the axes of first pc plot viewer should be \"AGE, WEIGHT, STARTED\"", () => axesShouldBe(page, el("first pc plot viewer"), "AGE, WEIGHT, STARTED"));
      await session.step(68, "And \"Title\" property of second pc plot viewer should not be \"Source Plot\"", () => propertyShouldNotBe(page, "Title", el("second pc plot viewer"), "Source Plot"));
      await session.step(69, "And \"Column Names\" property of second pc plot viewer should not be \"AGE, WEIGHT, STARTED\"", () => propertyShouldNotBe(page, "Column Names", el("second pc plot viewer"), "AGE, WEIGHT, STARTED"));
      await session.step(70, "When user picks \"Pick Up / Apply > Pick Up\" from the context menu of first pc plot viewer", () => pickFromContextMenu(page, "Pick Up / Apply > Pick Up", el("first pc plot viewer")));
      await session.step(71, "And user picks \"Pick Up / Apply > Apply\" from the context menu of second pc plot viewer", () => pickFromContextMenu(page, "Pick Up / Apply > Apply", el("second pc plot viewer")));
      await session.step(72, "Then properties of second pc plot viewer should be:", () => propertiesShouldBe(page, el("second pc plot viewer"), [["Title","Source Plot"],["Color","RACE"],["Legend Position","Left"],["Log Columns","AGE"],["Column Names","AGE, WEIGHT, STARTED"]]));
      await session.step(78, "And the axes of second pc plot viewer should be \"AGE, WEIGHT, STARTED\"", () => axesShouldBe(page, el("second pc plot viewer"), "AGE, WEIGHT, STARTED"));
      await session.step(79, "When user sets \"Column Names\" property of first pc plot viewer to \"AGE, HEIGHT, WEIGHT, STARTED\"", () => setProperty(page, "Column Names", el("first pc plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(80, "Then \"Column Names\" property of first pc plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => propertyShouldBe(page, "Column Names", el("first pc plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(81, "And \"Column Names\" property of second pc plot viewer should be \"AGE, WEIGHT, STARTED\"", () => propertyShouldBe(page, "Column Names", el("second pc plot viewer"), "AGE, WEIGHT, STARTED"));
      await session.step(82, "When user clicks on close icon of second pc plot viewer", () => clickOn(page, el("close icon of second pc plot viewer")));
      await session.step(83, "Then the open tableview should have 1 pc plot viewer", () => viewerCount(page, 1, "pc plot"));
      await session.step(84, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Column Names","AGE, HEIGHT, WEIGHT"],["Log Columns",""],["Legend Position","Auto"],["Title",""],["Show Title","false"],["Color",""]]));
      await session.step(91, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(92, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
