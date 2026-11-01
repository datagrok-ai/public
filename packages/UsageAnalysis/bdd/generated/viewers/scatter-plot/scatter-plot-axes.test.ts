/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-axes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, public#1882, public#2402, core#375]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll, passingAsRemembered, rememberPassing} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, legendLists, loadLayout, noBalloons, noErrors, painted, pickFromContextMenu, propertiesShouldBe, propertyShouldBe, readingAsRemembered, readingDoesNotRead, readingFinite, readingHigher, readingIs, rememberReading, repainted, saveLayoutToServer, setProperties, setProperty, showsFewerRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle, pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot axes, encodings and persistence", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-axes.feature", import.meta.url);
  test("Scatter plot axes, encodings and persistence", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot", "@realizes:public#1882", "@realizes:public#2402", "@realizes:core#375"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(17, "Then scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
    await run.scenario("The on-viewer selectors set the axes and the encodings", async () => {
      await session.step(20, "When user picks \"AGE\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "AGE", "x", el("scatter plot viewer")));
      await session.step(21, "And user picks \"HEIGHT\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "HEIGHT", "y", el("scatter plot viewer")));
      await session.step(22, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(23, "And user picks \"RACE\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "RACE", "color", el("scatter plot viewer")));
      await session.step(24, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(25, "And user picks \"WEIGHT\" in the \"size\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "size", el("scatter plot viewer")));
      await session.step(26, "And user sets \"Markers\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "SEX"));
      await session.step(27, "And user picks \"WEIGHT\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "WEIGHT", "x", el("scatter plot viewer")));
      await session.step(28, "And user picks \"AGE\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "AGE", "x", el("scatter plot viewer")));
      await session.step(29, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]);
      await session.step(35, "And X column input in scatter plot viewer should contain text \"AGE\"", () => shouldContainText(page, el("X column input in scatter plot viewer"), "AGE"));
      await session.step(36, "And Y column input in scatter plot viewer should contain text \"HEIGHT\"", () => shouldContainText(page, el("Y column input in scatter plot viewer"), "HEIGHT"));
      await session.step(37, "When user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
      await session.step(38, "Then Color column input in scatter plot viewer should contain text \"RACE\"", () => shouldContainText(page, el("Color column input in scatter plot viewer"), "RACE"));
      await session.step(39, "And Size column input in scatter plot viewer should contain text \"WEIGHT\"", () => shouldContainText(page, el("Size column input in scatter plot viewer"), "WEIGHT"));
      await session.step(40, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An explicit window, dropped by a new X column (public#1882), a logarithmic axis and an inverted one", async () => {
      await session.step(45, "When user sets \"X Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(46, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(48, "When user sets \"X Axis Type\" property of scatter plot viewer to \"linear\"", () => setProperty(page, "X Axis Type", el("scatter plot viewer"), "linear"));
      await session.step(49, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","20"],["xMax","60"]]), [["xMin","20"],["xMax","60"]]);
      await session.step(52, "Then the \"x axis min\" reading of scatter plot viewer should be 20", () => readingIs(page, "x axis min", el("scatter plot viewer"), 20));
      await session.step(53, "And the \"x axis max\" reading of scatter plot viewer should be 60", () => readingIs(page, "x axis max", el("scatter plot viewer"), 60));
      await session.step(54, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(55, "When user sets \"X\" property of scatter plot viewer to \"WEIGHT\"", () => setProperty(page, "X", el("scatter plot viewer"), "WEIGHT"));
      await session.step(56, "Then \"xMin\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "xMin", el("scatter plot viewer"), ""));
      await session.step(57, "And \"xMax\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "xMax", el("scatter plot viewer"), ""));
      await session.step(58, "And the \"x axis min\" reading of scatter plot viewer should not be \"20\"", () => readingDoesNotRead(page, "x axis min", el("scatter plot viewer"), "20"));
      await session.step(59, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","AGE"],["xMin","20"],["xMax","60"]]), [["X","AGE"],["xMin","20"],["xMax","60"]]);
      await session.step(63, "Then the \"x axis min\" reading of scatter plot viewer should be 20", () => readingIs(page, "x axis min", el("scatter plot viewer"), 20));
      await session.step(64, "When user sets \"Invert X Axis\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Invert X Axis", el("scatter plot viewer"), "true"));
      await session.step(65, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(66, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","60"],["xMax","20"]]), [["xMin","60"],["xMax","20"]]);
      await session.step(69, "Then the \"x axis min\" reading of scatter plot viewer should be 20", () => readingIs(page, "x axis min", el("scatter plot viewer"), 20));
      await session.step(70, "And the \"x axis max\" reading of scatter plot viewer should be 60", () => readingIs(page, "x axis max", el("scatter plot viewer"), 60));
      await session.step(71, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(72, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "When user sets \"X Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "X Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(75, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X Axis Type","logarithmic"],["Invert X Axis","true"],["xMin","60"],["xMax","20"]]), [["X Axis Type","logarithmic"],["Invert X Axis","true"],["xMin","60"],["xMax","20"]]);
      await session.step(80, "And scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
      await session.step(81, "And the \"x axis min\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis min", el("scatter plot viewer")));
      await session.step(82, "And the \"x axis max\" reading of scatter plot viewer should be a finite number", () => readingFinite(page, "x axis max", el("scatter plot viewer")));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
      await session.step(85, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin",""],["xMax",""],["Invert X Axis","false"],["X Axis Type","linear"]]), [["xMin",""],["xMax",""],["Invert X Axis","false"],["X Axis Type","linear"]]);
      await session.step(90, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["xMin",""],["xMax",""],["Invert X Axis","false"],["X Axis Type","linear"]]), [["xMin",""],["xMax",""],["Invert X Axis","false"],["X Axis Type","linear"]]);
      await session.step(95, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A datetime X axis takes the axis type away and offers the time unit instead", async () => {
      await session.step(99, "When user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
      await session.step(100, "And user sets \"X\" property of scatter plot viewer to \"STARTED\"", () => setProperty(page, "X", el("scatter plot viewer"), "STARTED"));
      await session.step(101, "Then \"X Axis Type\" property should be disabled", () => shouldBe(page, el("\"X Axis Type\" property"), "disabled"));
      await session.step(102, "And \"X Map\" property should be enabled", () => shouldBe(page, el("\"X Map\" property"), "enabled"));
      await session.step(103, "When user sets \"X\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "X", el("scatter plot viewer"), "AGE"));
      await session.step(104, "Then \"X Axis Type\" property should be enabled", () => shouldBe(page, el("\"X Axis Type\" property"), "enabled"));
      await session.step(105, "And \"X Map\" property should be disabled", () => shouldBe(page, el("\"X Map\" property"), "disabled"));
      await session.step(106, "And \"X\" property of scatter plot viewer should be \"AGE\"", () => propertyShouldBe(page, "X", el("scatter plot viewer"), "AGE"));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("One column serves both axes", async () => {
      await session.step(110, "When user sets \"Y\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Y", el("scatter plot viewer"), "AGE"));
      await session.step(111, "Then properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","AGE"]]), [["X","AGE"],["Y","AGE"]]);
      await session.step(114, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(115, "And scatter plot viewer should be painted", () => painted(page, el("scatter plot viewer")));
      await session.step(116, "When user sets \"Y\" property of scatter plot viewer to \"HEIGHT\"", () => setProperty(page, "Y", el("scatter plot viewer"), "HEIGHT"));
      await session.step(117, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(118, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout round-trip on the server brings the five columns back, and an explicit X window (public#2402)", async () => {
      await session.step(122, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin","30"],["xMax","60"]]), [["xMin","30"],["xMax","60"]]);
      await session.step(125, "Then the \"x axis min\" reading of scatter plot viewer should be 30", () => readingIs(page, "x axis min", el("scatter plot viewer"), 30));
      await session.step(126, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(127, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(128, "And user sets \"Color\" property of scatter plot viewer to \"\"", () => setProperty(page, "Color", el("scatter plot viewer"), ""));
      await session.step(129, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(130, "When user loads the saved layout", () => loadLayout(page));
      await session.step(131, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(132, "And histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(133, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]);
      await session.step(139, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(140, "And the \"x axis min\" reading of scatter plot viewer should be 30", () => readingIs(page, "x axis min", el("scatter plot viewer"), 30));
      await session.step(141, "And the \"x axis max\" reading of scatter plot viewer should be 60", () => readingIs(page, "x axis max", el("scatter plot viewer"), 60));
      await session.step(142, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["xMin",""],["xMax",""]]), [["xMin",""],["xMax",""]]);
      await session.step(145, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A project round-trip brings them back too, with the window a zoom left (core#375)", async () => {
      await session.step(150, "When user drags the min handle of the \"x\" range slider of scatter plot viewer by 40 pixels", () => dragRangeHandle(page, "min", "x", el("scatter plot viewer"), 40));
      await session.step(151, "Then the \"x axis min\" reading of scatter plot viewer should be higher than before", () => readingHigher(page, "x axis min", el("scatter plot viewer")));
      await session.step(152, "When user remembers the \"x axis min\" reading of scatter plot viewer", () => rememberReading(page, "x axis min", el("scatter plot viewer")));
      await session.step(153, "And user remembers how many rows pass the filter", () => rememberPassing(page));
      await session.step(154, "And user saves the current view as project \"zz-scatter-plot-axes\"", () => saveAsProject(page, "zz-scatter-plot-axes"));
      await session.step(155, "And user closes all views", () => closeAllViews(page));
      await session.step(156, "And user opens the \"zz-scatter-plot-axes\" project", () => openProject(page, "zz-scatter-plot-axes"));
      await session.step(157, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(158, "And properties of scatter plot viewer should be:", () => propertiesShouldBe(page, el("scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]), [["X","AGE"],["Y","HEIGHT"],["Color","RACE"],["Size","WEIGHT"],["Markers","SEX"]]);
      await session.step(164, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(165, "And the \"x axis min\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis min", el("scatter plot viewer")));
      await session.step(166, "And as many rows as remembered should pass the filter", () => passingAsRemembered(page));
      await session.step(167, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(168, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(169, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(170, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
