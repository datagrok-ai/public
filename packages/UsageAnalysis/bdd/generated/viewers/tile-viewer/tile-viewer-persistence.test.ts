/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addTileViewerWith, deleteValueField} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, loadLayout, noErrors, propertiesShouldBe, readingAsRemembered, readingHigher, readingIs, readingReads, rememberReading, saveLayout, saveLayoutToServer, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu, readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer persistence", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-persistence.feature", import.meta.url);
  test("Tile viewer persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.tile-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"],["Lanes","Black, Asian"],["Show Title","true"],["Title","Patient cards"]]));
    await session.step(28, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
    await session.step(29, "And the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
    await session.step(30, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
    await run.scenario("The configured viewer is the one on screen", async () => {
      await session.step(33, "Then title of tile viewer should have text \"Patient cards\"", () => shouldHaveText(page, el("title of tile viewer"), "Patient cards"));
      await session.step(34, "And the \"lanes list\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lanes list", el("tile viewer"), "Black, Asian"));
      await session.step(35, "And tile viewer should have a \"lane header Black\" area", () => hasArea(page, el("tile viewer"), "lane header Black"));
      await session.step(36, "And tile viewer should not have a \"lane Caucasian\" area", () => hasNoArea(page, el("tile viewer"), "lane Caucasian"));
      await session.step(37, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(38, "And the \"lane of row 10\" reading of tile viewer should be \"Asian\"", () => readingReads(page, "lane of row 10", el("tile viewer"), "Asian"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A scrolled lane keeps its position when another viewer is docked beside it", async () => {
      await session.step(42, "Then the \"scroll of lane Black\" reading of tile viewer should be 0", () => readingIs(page, "scroll of lane Black", el("tile viewer"), 0));
      await session.step(43, "When user scrolls the mouse wheel down over the \"lane content Black\" area of tile viewer", () => wheelOverArea(page, "down", "lane content Black", el("tile viewer")));
      await session.step(44, "Then the \"scroll of lane Black\" reading of tile viewer should be higher than before", () => readingHigher(page, "scroll of lane Black", el("tile viewer")));
      await session.step(45, "When user remembers the \"scroll of lane Black\" reading of tile viewer", () => rememberReading(page, "scroll of lane Black", el("tile viewer")));
      await session.step(46, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(47, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
      await session.step(48, "And the \"scroll of lane Black\" reading of tile viewer should be as remembered", () => readingAsRemembered(page, "scroll of lane Black", el("tile viewer")));
      await session.step(49, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(50, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved on the server brings the lanes and the viewer set back", async () => {
      await session.step(54, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(55, "And user clicks on close icon of tile viewer", () => clickOn(page, el("close icon of tile viewer")));
      await session.step(56, "Then tile viewer should be absent", () => shouldBe(page, el("tile viewer"), "absent"));
      await session.step(57, "When user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(58, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(59, "When user loads the saved layout", () => loadLayout(page));
      await session.step(60, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
      await session.step(61, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(62, "And properties of tile viewer should be:", () => propertiesShouldBe(page, el("tile viewer"), [["Lanes Column Name","RACE"],["Lanes","Black, Asian"],["Title","Patient cards"]]));
      await session.step(66, "And the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(67, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(68, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A designed field set survives a layout round-trip", async () => {
      await session.step(72, "When user picks \"Edit Form...\" from the viewer menu of tile viewer", () => pickFromViewerMenu(page, "Edit Form...", el("tile viewer")));
      await session.step(73, "Then form designer should be visible", () => shouldBe(page, el("form designer"), "visible"));
      await session.step(74, "When user deletes the \"WEIGHT\" value field in the form designer", () => deleteValueField(page, "WEIGHT"));
      await session.step(75, "And user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
      await session.step(76, "Then form designer should be absent", () => shouldBe(page, el("form designer"), "absent"));
      await session.step(77, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(78, "And the \"fields\" reading of tile viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "fields", el("tile viewer"), "WEIGHT"));
      await session.step(79, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(80, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(81, "And user clicks on close icon of tile viewer", () => clickOn(page, el("close icon of tile viewer")));
      await session.step(82, "Then tile viewer should be absent", () => shouldBe(page, el("tile viewer"), "absent"));
      await session.step(83, "When user loads the saved layout", () => loadLayout(page));
      await session.step(84, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
      await session.step(85, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(86, "And the \"fields\" reading of tile viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "fields", el("tile viewer"), "WEIGHT"));
      await session.step(87, "And the \"fields\" reading of tile viewer should contain \"AGE\"", () => readingContains(page, "fields", el("tile viewer"), "AGE"));
      await session.step(88, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(89, "And the \"auto generate\" reading of tile viewer should be \"false\"", () => readingReads(page, "auto generate", el("tile viewer"), "false"));
      await session.step(90, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project saved, closed and reopened brings all of it back", async () => {
      await session.step(94, "When user saves the current view as project \"bdd tile viewer round trip\"", () => saveAsProject(page, "bdd tile viewer round trip"));
      await session.step(95, "And user closes all views", () => closeAllViews(page));
      await session.step(96, "And user opens the \"bdd tile viewer round trip\" project", () => openProject(page, "bdd tile viewer round trip"));
      await session.step(97, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
      await session.step(98, "And properties of tile viewer should be:", () => propertiesShouldBe(page, el("tile viewer"), [["Lanes Column Name","RACE"],["Lanes","Black, Asian"],["Title","Patient cards"]]));
      await session.step(102, "And the \"lanes\" reading of tile viewer should be 2", () => readingIs(page, "lanes", el("tile viewer"), 2));
      await session.step(103, "And the \"lane names\" reading of tile viewer should be \"Black, Asian\"", () => readingReads(page, "lane names", el("tile viewer"), "Black, Asian"));
      await session.step(104, "And the \"fields shown\" reading of tile viewer should be 9", () => readingIs(page, "fields shown", el("tile viewer"), 9));
      await session.step(105, "And the \"fields\" reading of tile viewer should not contain \"WEIGHT\"", () => readingNotContains(page, "fields", el("tile viewer"), "WEIGHT"));
      await session.step(106, "And the \"form designed\" reading of tile viewer should be \"true\"", () => readingReads(page, "form designed", el("tile viewer"), "true"));
      await session.step(107, "And the \"lane of row 101\" reading of tile viewer should be \"Black\"", () => readingReads(page, "lane of row 101", el("tile viewer"), "Black"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
