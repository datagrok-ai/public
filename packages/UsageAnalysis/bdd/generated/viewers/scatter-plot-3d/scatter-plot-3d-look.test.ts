/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot-3d/scatter-plot-3d-look.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot3d]
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
import {clickOn, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaNarrower, areaTaller, areaWider, errorBalloonText, noBalloons, noErrors, pickColorSwatch, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingDiffers, readingIs, readingNotAsRemembered, rememberReading, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("3D scatter plot look", () => {
  const session = feature(test, "features/viewers/scatter-plot-3d/scatter-plot-3d-look.feature", import.meta.url);
  test("3D scatter plot look", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
    await session.step(22, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await run.scenario("Y and Z axes go logarithmic and come back", async () => {
      await session.step(25, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(26, "And user sets \"Y Axis Type\" property of 3d scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("3d scatter plot viewer"), "logarithmic"));
      await session.step(27, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(28, "When user sets \"Z Axis Type\" property of 3d scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Z Axis Type", el("3d scatter plot viewer"), "logarithmic"));
      await session.step(29, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(30, "And the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(32, "When user sets \"Y Axis Type\" property of 3d scatter plot viewer to \"linear\"", () => setProperty(page, "Y Axis Type", el("3d scatter plot viewer"), "linear"));
      await session.step(33, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(34, "When user sets \"Z Axis Type\" property of 3d scatter plot viewer to \"linear\"", () => setProperty(page, "Z Axis Type", el("3d scatter plot viewer"), "linear"));
      await session.step(35, "Then the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Size by a number scales the markers and Marker Size Scaling redraws them", async () => {
      await session.step(39, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(40, "And user sets \"Size\" property of 3d scatter plot viewer to \"WEIGHT\"", () => setProperty(page, "Size", el("3d scatter plot viewer"), "WEIGHT"));
      await session.step(41, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(42, "When user sets \"Marker Size Scaling\" property of 3d scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Marker Size Scaling", el("3d scatter plot viewer"), "logarithmic"));
      await session.step(43, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(44, "When user sets \"Marker Size Scaling\" property of 3d scatter plot viewer to \"linear\"", () => setProperty(page, "Marker Size Scaling", el("3d scatter plot viewer"), "linear"));
      await session.step(45, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(46, "When user sets \"Size\" property of 3d scatter plot viewer to \"\"", () => setProperty(page, "Size", el("3d scatter plot viewer"), ""));
      await session.step(47, "Then the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A label column draws labels and clearing it takes them away", async () => {
      await session.step(51, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(52, "And user sets \"Label\" property of 3d scatter plot viewer to \"SEX\"", () => setProperty(page, "Label", el("3d scatter plot viewer"), "SEX"));
      await session.step(53, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(54, "When user sets \"Label\" property of 3d scatter plot viewer to \"\"", () => setProperty(page, "Label", el("3d scatter plot viewer"), ""));
      await session.step(55, "Then the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Random marker rotation turns the markers and switching it off puts them back", async () => {
      await session.step(59, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(60, "And user sets \"Marker Random Rotation\" property of 3d scatter plot viewer to \"true\"", () => setProperty(page, "Marker Random Rotation", el("3d scatter plot viewer"), "true"));
      await session.step(61, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(62, "When user sets \"Marker Random Rotation\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Marker Random Rotation", el("3d scatter plot viewer"), "false"));
      await session.step(63, "Then the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Back Color picked in the context panel repaints the scene behind the markers", async () => {
      await session.step(67, "When user remembers the \"scene signature\" reading of 3d scatter plot viewer", () => rememberReading(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(68, "And user picks \"Properties...\" from the context menu of 3d scatter plot viewer", () => pickFromContextMenu(page, "Properties...", el("3d scatter plot viewer")));
      await session.step(69, "Then the context panel should show \"3d scatter plot\"", () => contextPanelShows(page, "3d scatter plot"));
      await session.step(70, "Given \"Style\" category in context panel is expanded", () => isExpanded(page, el("\"Style\" category in context panel")));
      await session.step(71, "When user clicks on editor of \"Back Color\" property in context panel", () => clickOn(page, el("editor of \"Back Color\" property in context panel")));
      await session.step(72, "And user picks the color \"#FFA500\" in the color picker", () => pickColorSwatch(page, "#FFA500"));
      await session.step(73, "Then \"Back Color\" property of 3d scatter plot viewer should be \"#FFA500\"", () => propertyShouldBe(page, "Back Color", el("3d scatter plot viewer"), "#FFA500"));
      await session.step(74, "And the \"scene signature\" reading of 3d scatter plot viewer should not be as remembered", () => readingNotAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(75, "When user sets \"Back Color\" property of 3d scatter plot viewer to \"#FFFFFF\"", () => setProperty(page, "Back Color", el("3d scatter plot viewer"), "#FFFFFF"));
      await session.step(76, "Then the \"scene signature\" reading of 3d scatter plot viewer should be as remembered", () => readingAsRemembered(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Hidden column selectors leave the plot and give their strip to the scene", async () => {
      await session.step(80, "When user sets \"Show X Selector\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show X Selector", el("3d scatter plot viewer"), "false"));
      await session.step(81, "And user sets \"Show X Selector\" property of 3d scatter plot viewer to \"true\"", () => setProperty(page, "Show X Selector", el("3d scatter plot viewer"), "true"));
      await session.step(82, "Then X column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("X column input in 3d scatter plot viewer"), "visible"));
      await session.step(83, "When user sets \"Show X Selector\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show X Selector", el("3d scatter plot viewer"), "false"));
      await session.step(84, "Then X column input in 3d scatter plot viewer should be hidden", () => shouldBe(page, el("X column input in 3d scatter plot viewer"), "hidden"));
      await session.step(85, "And the \"view\" area of 3d scatter plot viewer should be taller than before", () => areaTaller(page, "view", el("3d scatter plot viewer")));
      await session.step(86, "When user sets \"Show Color Selector\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Color Selector", el("3d scatter plot viewer"), "false"));
      await session.step(87, "Then Color column input in 3d scatter plot viewer should be hidden", () => shouldBe(page, el("Color column input in 3d scatter plot viewer"), "hidden"));
      await session.step(88, "And the \"view\" area of 3d scatter plot viewer should be taller than before", () => areaTaller(page, "view", el("3d scatter plot viewer")));
      await session.step(89, "When user sets \"Show Y Selector\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Y Selector", el("3d scatter plot viewer"), "false"));
      await session.step(90, "Then Y column input in 3d scatter plot viewer should be hidden", () => shouldBe(page, el("Y column input in 3d scatter plot viewer"), "hidden"));
      await session.step(91, "And the \"view\" area of 3d scatter plot viewer should be wider than before", () => areaWider(page, "view", el("3d scatter plot viewer")));
      await session.step(92, "When user sets \"Show Z Selector\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Z Selector", el("3d scatter plot viewer"), "false"));
      await session.step(93, "Then Z column input in 3d scatter plot viewer should be hidden", () => shouldBe(page, el("Z column input in 3d scatter plot viewer"), "hidden"));
      await session.step(94, "And the \"view\" area of 3d scatter plot viewer should be wider than before", () => areaWider(page, "view", el("3d scatter plot viewer")));
      await session.step(95, "When user sets properties of 3d scatter plot viewer:", () => setProperties(page, el("3d scatter plot viewer"), [["Show X Selector","true"],["Show Y Selector","true"],["Show Z Selector","true"],["Show Color Selector","true"]]), [["Show X Selector","true"],["Show Y Selector","true"],["Show Z Selector","true"],["Show Color Selector","true"]]);
      await session.step(100, "Then X column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("X column input in 3d scatter plot viewer"), "visible"));
      await session.step(101, "And Y column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Y column input in 3d scatter plot viewer"), "visible"));
      await session.step(102, "And Z column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Z column input in 3d scatter plot viewer"), "visible"));
      await session.step(103, "And Color column input in 3d scatter plot viewer should be visible", () => shouldBe(page, el("Color column input in 3d scatter plot viewer"), "visible"));
      await session.step(104, "And the \"view\" area of 3d scatter plot viewer should be narrower than before", () => areaNarrower(page, "view", el("3d scatter plot viewer")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis over zero values says why it cannot be drawn", async () => {
      await session.step(108, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(109, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
      await session.step(110, "Then Y column input in 3d scatter plot viewer should contain text \"turbo\"", () => shouldContainText(page, el("Y column input in 3d scatter plot viewer"), "turbo"));
      await session.step(111, "When user sets \"Y Axis Type\" property of 3d scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("3d scatter plot viewer"), "logarithmic"));
      await session.step(112, "Then an error balloon containing \"'turbo' contains negative or zero values and can't be visualized in the log scale\" should have been shown", () => errorBalloonText(page, "'turbo' contains negative or zero values and can't be visualized in the log scale"));
      await session.step(113, "When user sets \"Y Axis Type\" property of 3d scatter plot viewer to \"linear\"", () => setProperty(page, "Y Axis Type", el("3d scatter plot viewer"), "linear"));
      await session.step(114, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 30", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 30));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
