/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-persistence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {slicesOffCentre} from '../../../bindings/pie-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {categoricalColorIs, colorCategorical, colorCodedCategorically, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, areaColor, areaNotColor, legendItemColor, loadLayout, noErrors, propertiesShouldBe, readingBetween, readingIs, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart persistence", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-persistence.feature", import.meta.url);
  test("Pie chart persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"],["Segment Angle Aggr Type","sum"],["Show Value","true"],["Start Angle","45"],["Shift","5"],["Legend Visibility","Always"],["Show Title","true"],["Title","Race mix"]]));
    await session.step(25, "And user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Asian","#9467BD"]]));
    await session.step(27, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(28, "And the \"start angle of Asian\" reading of pie chart viewer should be 45", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 45));
    await session.step(29, "And the \"share of Caucasian\" reading of pie chart viewer should be between 89.55 and 89.6", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 89.55, 89.6));
    await run.scenario("The configured chart is the one on screen", async () => {
      await session.step(32, "Then title of pie chart viewer should have text \"Race mix\"", () => shouldHaveText(page, el("title of pie chart viewer"), "Race mix"));
      await session.step(33, "And the slices of pie chart viewer should sit 5 pixels off the centre", () => slicesOffCentre(page, el("pie chart viewer"), 5));
      await session.step(34, "And \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(35, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#9467BD"));
      await session.step(36, "And the \"pie\" area of pie chart viewer should contain the color \"#9467BD\"", () => areaColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved on the server restores the configuration and the viewer set", async () => {
      await session.step(40, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(41, "And user clicks on close icon of pie chart viewer", () => clickOn(page, el("close icon of pie chart viewer")));
      await session.step(42, "Then pie chart viewer should be absent", () => shouldBe(page, el("pie chart viewer"), "absent"));
      await session.step(43, "When user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(44, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(45, "When user loads the saved layout", () => loadLayout(page));
      await session.step(46, "Then pie chart viewer should be visible", () => shouldBe(page, el("pie chart viewer"), "visible"));
      await session.step(47, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(48, "And properties of pie chart viewer should be:", () => propertiesShouldBe(page, el("pie chart viewer"), [["Category","RACE"],["Segment Angle Aggr Type","sum"],["Show Value","true"],["Start Angle","45"],["Shift","5"],["Title","Race mix"]]));
      await session.step(55, "And the \"start angle of Asian\" reading of pie chart viewer should be 45", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 45));
      await session.step(56, "And the \"share of Caucasian\" reading of pie chart viewer should be between 89.55 and 89.6", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 89.55, 89.6));
      await session.step(57, "And the slices of pie chart viewer should sit 5 pixels off the centre", () => slicesOffCentre(page, el("pie chart viewer"), 5));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The custom category colour comes back with the layout", async () => {
      await session.step(61, "Then \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(62, "And the categorical color of \"Asian\" in \"RACE\" column should be \"#9467BD\"", () => categoricalColorIs(page, "Asian", "RACE", "#9467BD"));
      await session.step(63, "And the \"pie\" area of pie chart viewer should contain the color \"#9467BD\"", () => areaColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project saved, closed and reopened brings all of it back", async () => {
      await session.step(67, "When user saves the current view as project \"bdd pie chart round trip\"", () => saveAsProject(page, "bdd pie chart round trip"));
      await session.step(68, "And user closes all views", () => closeAllViews(page));
      await session.step(69, "And user opens the \"bdd pie chart round trip\" project", () => openProject(page, "bdd pie chart round trip"));
      await session.step(70, "Then pie chart viewer should be visible", () => shouldBe(page, el("pie chart viewer"), "visible"));
      await session.step(71, "And properties of pie chart viewer should be:", () => propertiesShouldBe(page, el("pie chart viewer"), [["Category","RACE"],["Segment Angle Aggr Type","sum"],["Show Value","true"],["Start Angle","45"],["Shift","5"],["Title","Race mix"]]));
      await session.step(78, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(79, "And the \"start angle of Asian\" reading of pie chart viewer should be 45", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 45));
      await session.step(80, "And the \"share of Caucasian\" reading of pie chart viewer should be between 89.55 and 89.6", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 89.55, 89.6));
      await session.step(81, "And \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(82, "And the \"pie\" area of pie chart viewer should contain the color \"#9467BD\"", () => areaColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The coloring goes back where it was found", async () => {
      await session.step(86, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Asian","#1F77B4"]]));
      await session.step(88, "Then the \"pie\" area of pie chart viewer should contain the color \"#1F77B4\"", () => areaColor(page, "pie", el("pie chart viewer"), "#1F77B4"));
      await session.step(89, "And the \"pie\" area of pie chart viewer should not contain the color \"#9467BD\"", () => areaNotColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(90, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(91, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
