/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-settings-ladder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {zoomValueAxis} from '../../../bindings/box-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, loadLayout, narrowerRange, noErrors, notRepainted, propertiesShouldBe, propertyShouldBe, rangeWithinColumn, rememberRange, rememberedRange, repainted, sameRange, saveLayout, setProperties, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot settings ladder", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-settings-ladder.feature", import.meta.url);
  test("Box plot settings ladder", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"]]));
    await run.scenario("A datetime value disables Axis Type", async () => {
      await session.step(16, "When user sets \"Value\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Value", el("box plot viewer"), "STARTED"));
      await session.step(17, "And user clicks on settings icon of box plot viewer", () => clickOn(page, el("settings icon of box plot viewer")));
      await session.step(18, "Then \"Axis Type\" property in context panel should be disabled", () => shouldBe(page, el("\"Axis Type\" property in context panel"), "disabled"));
      await session.step(19, "When user sets \"Value\" property of box plot viewer to \"AGE\"", () => setProperty(page, "Value", el("box plot viewer"), "AGE"));
      await session.step(20, "Then \"Axis Type\" property in context panel should be enabled", () => shouldBe(page, el("\"Axis Type\" property in context panel"), "enabled"));
    });
    await run.scenario("Category 1 sets the marker color", async () => {
      await session.step(23, "Then \"Category 1\" property of box plot viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "DIS_POP"));
      await session.step(24, "And \"Marker Color Column\" property of box plot viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "DIS_POP"));
      await session.step(25, "When user sets \"Category 1\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 1", el("box plot viewer"), "SEX"));
      await session.step(26, "Then \"Category 1\" property of box plot viewer should be \"SEX\"", () => propertyShouldBe(page, "Category 1", el("box plot viewer"), "SEX"));
      await session.step(27, "And \"Marker Color Column\" property of box plot viewer should be \"SEX\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "SEX"));
    });
    await run.scenario("An explicit coloring survives a value change", async () => {
      await session.step(30, "When user sets \"Marker Color Column\" property of box plot viewer to \"HEIGHT\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "HEIGHT"));
      await session.step(31, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Invert Color Scheme","true"],["Color Min","20"],["Color Max","80"]]));
      await session.step(35, "And user sets \"Category 2\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 2", el("box plot viewer"), "RACE"));
      await session.step(36, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Minor Categories","true"],["Show All Categories","true"]]));
      await session.step(39, "And user sets \"Value\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Value", el("box plot viewer"), "WEIGHT"));
      await session.step(40, "And user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(41, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(42, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Value","WEIGHT"],["Category 1","SEX"],["Category 2","RACE"],["Marker Color Column","HEIGHT"],["Invert Color Scheme","true"],["Color Min","20"],["Color Max","80"],["Show Minor Categories","true"],["Show All Categories","true"]]));
    });
    await run.scenario("Value limits and the log axis", async () => {
      await session.step(54, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value Min","20"],["Value Max","60"]]));
      await session.step(57, "Then properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Value Min","20"],["Value Max","60"]]));
      await session.step(60, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value Min",""],["Value Max",""]]));
      await session.step(63, "And user sets \"Axis Type\" property of box plot viewer to \"logarithmic\"", () => setProperty(page, "Axis Type", el("box plot viewer"), "logarithmic"));
      await session.step(64, "Then \"Axis Type\" property of box plot viewer should be \"logarithmic\"", () => propertyShouldBe(page, "Axis Type", el("box plot viewer"), "logarithmic"));
      await session.step(65, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(66, "And no errors should have been logged", () => noErrors(page));
      await session.step(67, "And the value range of box plot viewer should lie within \"WEIGHT\" column", () => rangeWithinColumn(page, el("box plot viewer"), "WEIGHT"));
      await session.step(68, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Invert Y Axis","true"],["Plot Style","violin"]]));
      await session.step(71, "Then properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Invert Y Axis","true"],["Plot Style","violin"]]));
      await session.step(74, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
    });
    await run.scenario("A zoom survives a coloring change", async () => {
      await session.step(77, "When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await session.step(78, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(79, "When user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(80, "Then \"Marker Color Column\" property of box plot viewer should be \"SEX\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(81, "And box plot viewer should show the same value range as before", () => sameRange(page, el("box plot viewer")));
    });
    await run.scenario("Group comparison with a control and a covariate", async () => {
      await session.step(84, "When user sets \"Show Group Comparison\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "true"));
      await session.step(85, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Control Comparisons","true"],["Control Group","F"]]));
      await session.step(88, "Then \"Control Group\" property of box plot viewer should be \"F\"", () => propertyShouldBe(page, "Control Group", el("box plot viewer"), "F"));
      await session.step(89, "When user sets \"Adjust By\" property of box plot viewer to \"HEIGHT\"", () => setProperty(page, "Adjust By", el("box plot viewer"), "HEIGHT"));
      await session.step(90, "Then \"Adjust By\" property of box plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "Adjust By", el("box plot viewer"), "HEIGHT"));
      await session.step(91, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
    });
    await run.scenario("The ladder survives a layout round-trip", async () => {
      await session.step(94, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Adjust By",""],["Control Comparisons","false"],["Show Group Comparison","false"]]));
      await session.step(98, "And user saves the layout of the current table view", () => saveLayout(page));
      await session.step(99, "And user clicks on close icon of box plot viewer", () => clickOn(page, el("close icon of box plot viewer")));
      await session.step(100, "Then box plot viewer should be absent", () => shouldBe(page, el("box plot viewer"), "absent"));
      await session.step(101, "When user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
      await session.step(102, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(103, "When user loads the saved layout", () => loadLayout(page));
      await session.step(104, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
      await session.step(105, "And scatter plot viewer should be absent", () => shouldBe(page, el("scatter plot viewer"), "absent"));
      await session.step(106, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Value","WEIGHT"],["Category 1","SEX"],["Category 2","RACE"],["Show Minor Categories","true"],["Show All Categories","true"],["Marker Color Column","SEX"],["Invert Color Scheme","true"],["Color Min","20"],["Color Max","80"],["Axis Type","logarithmic"],["Invert Y Axis","true"],["Plot Style","violin"]]));
    });
    await run.scenario("The ladder survives a project round-trip", async () => {
      await session.step(121, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Group Comparison","true"],["Control Comparisons","true"],["Control Group","F"],["Adjust By","HEIGHT"]]));
      await session.step(126, "And user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
      await session.step(127, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
      await session.step(128, "When user remembers the value range of box plot viewer", () => rememberRange(page, el("box plot viewer")));
      await session.step(129, "And user saves the current view as project \"bdd box plot ladder\"", () => saveAsProject(page, "bdd box plot ladder"));
      await session.step(130, "And user closes all views", () => closeAllViews(page));
      await session.step(131, "And user opens the \"bdd box plot ladder\" project", () => openProject(page, "bdd box plot ladder"));
      await session.step(132, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
      await session.step(133, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Value","WEIGHT"],["Category 1","SEX"],["Category 2","RACE"],["Marker Color Column","SEX"],["Invert Color Scheme","true"],["Color Min","20"],["Color Max","80"],["Axis Type","logarithmic"],["Invert Y Axis","true"],["Plot Style","violin"],["Show Group Comparison","true"],["Control Group","F"],["Adjust By","HEIGHT"]]));
      await session.step(147, "And box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
    });
    run.finish();
  });
});
