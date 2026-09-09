/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-axes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {pickColumn} from '../../../bindings/scatter-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, legendLists, loadLayout, noBalloons, noErrors, painted, propertiesShouldBe, propertyShouldBe, readingIs, repainted, saveLayoutToServer, setProperties, setProperty, showsFewerRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot axes, encodings and persistence", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-axes.feature", import.meta.url);
  test("Scatter plot axes, encodings and persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(14, "Then scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await run.scenario("The on-viewer selectors set the axes and the encodings", async () => {
      await session.step(17, "When user picks \"AGE\" in the X column selector of scatter plot viewer", () => pickColumn(page, "AGE", "X"));
      await session.step(18, "And user picks \"HEIGHT\" in the Y column selector of scatter plot viewer", () => pickColumn(page, "HEIGHT", "Y"));
      await session.step(19, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(20, "And user picks \"RACE\" in the Color column selector of scatter plot viewer", () => pickColumn(page, "RACE", "Color"));
      await session.step(21, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(22, "And user picks \"WEIGHT\" in the Size column selector of scatter plot viewer", () => pickColumn(page, "WEIGHT", "Size"));
      await session.step(23, "And user sets \"Markers\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "SEX"));
      await session.step(24, "And user picks \"WEIGHT\" in the X column selector of scatter plot viewer", () => pickColumn(page, "WEIGHT", "X"));
      await session.step(25, "And user picks \"AGE\" in the X column selector of scatter plot viewer", () => pickColumn(page, "AGE", "X"));
      await session.step(26, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]));
      await session.step(32, "And X column input in scatter plot viewer should contain text \"AGE\"", () => shouldContainText(page, el("X column input in scatter plot viewer"), "AGE"));
      await session.step(33, "And Y column input in scatter plot viewer should contain text \"HEIGHT\"", () => shouldContainText(page, el("Y column input in scatter plot viewer"), "HEIGHT"));
      await session.step(34, "When user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(35, "Then Color column input in scatter plot viewer should contain text \"RACE\"", () => shouldContainText(page, el("Color column input in scatter plot viewer"), "RACE"));
      await session.step(36, "And Size column input in scatter plot viewer should contain text \"WEIGHT\"", () => shouldContainText(page, el("Size column input in scatter plot viewer"), "WEIGHT"));
      await session.step(37, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An explicit window, a logarithmic axis and an inverted one", async () => {
      await session.step(41, "When user sets \"X Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(42, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(44, "When user sets \"X Axis Type\" property of scatter plot viewer to \"linear\"", () => setProperty(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(45, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","20"],["xMax","60"]]));
      await session.step(48, "Then the \"x axis min\" reading of scatter plot viewer should be 20", () => readingIs(page, "x axis min", el("scatter plot viewer"), 20));
      await session.step(49, "And the \"x axis max\" reading of scatter plot viewer should be 60", () => readingIs(page, "x axis max", el("scatter plot viewer"), 60));
      await session.step(50, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(51, "When user sets \"Invert X Axis\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Invert X Axis", el("scatter plot viewer"), "true"));
      await session.step(52, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(53, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","60"],["xMax","20"]]));
      await session.step(56, "Then the \"x axis min\" reading of scatter plot viewer should be 20", () => readingIs(page, "x axis min", el("scatter plot viewer"), 20));
      await session.step(57, "And the \"x axis max\" reading of scatter plot viewer should be 60", () => readingIs(page, "x axis max", el("scatter plot viewer"), 60));
      await session.step(58, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(59, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin",""],["xMax",""],["Invert X Axis","false"]]));
      await session.step(65, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["xMin",""],["xMax",""],["Invert X Axis","false"],["X Axis Type","linear"]]));
      await session.step(70, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A datetime X axis takes the axis type away and offers the time unit instead", async () => {
      await session.step(74, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(75, "And user sets \"X\" property of scatter plot viewer to \"STARTED\"", () => setProperty(page, "X", el("scatter plot viewer"), "STARTED"));
      await session.step(76, "Then \"X Axis Type\" property should be disabled", () => shouldBe(page, el("\"X Axis Type\" property"), "disabled"));
      await session.step(77, "And \"X Map\" property should be enabled", () => shouldBe(page, el("\"X Map\" property"), "enabled"));
      await session.step(78, "When user sets \"X\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "X", el("scatter plot viewer"), "AGE"));
      await session.step(79, "Then \"X Axis Type\" property should be enabled", () => shouldBe(page, el("\"X Axis Type\" property"), "enabled"));
      await session.step(80, "And \"X Map\" property should be disabled", () => shouldBe(page, el("\"X Map\" property"), "disabled"));
      await session.step(81, "And \"X\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "AGE"));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One column serves both axes", async () => {
      await session.step(85, "When user sets \"Y\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Y", el("scatter plot viewer"), "AGE"));
      await session.step(86, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","AGE"]]));
      await session.step(89, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(90, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(91, "When user sets \"Y\" property of scatter plot viewer to \"HEIGHT\"", () => setProperty(page, "Y", el("scatter plot viewer"), "HEIGHT"));
      await session.step(92, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout round-trip on the server brings the five columns back", async () => {
      await session.step(96, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(97, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(98, "And user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(99, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(100, "When user loads the saved layout", () => loadLayout(page));
      await session.step(101, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(102, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(103, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]));
      await session.step(109, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings them back too", async () => {
      await session.step(113, "When user saves the current view as project \"zz-scatter-plot-axes\"", () => saveAsProject(page, "zz-scatter-plot-axes"));
      await session.step(114, "And user closes all views", () => closeAllViews(page));
      await session.step(115, "And user opens the \"zz-scatter-plot-axes\" project", () => openProject(page, "zz-scatter-plot-axes"));
      await session.step(116, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(117, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]));
      await session.step(123, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(124, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
