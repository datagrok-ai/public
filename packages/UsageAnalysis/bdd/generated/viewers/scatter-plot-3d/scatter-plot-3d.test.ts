/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot-3d/scatter-plot-3d.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot3d]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldHaveItems} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, filterPasses, filterTo, hasCurrentRow, resetFilter, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, clickAreaHolding, dragAcrossArea, hasArea, hoverArea, legendSide, noBalloons, noErrors, pickFromContextMenu, pointerAway, propertiesShouldBe, propertyShouldBe, readingDiffers, readingHigher, readingIs, readingLower, setProperties, setProperty, takeSnapshot, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("3D scatter plot", () => {
  const session = feature(test, "features/viewers/scatter-plot-3d/scatter-plot-3d.feature", import.meta.url);
  test("3D scatter plot", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot3d"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a 3d scatter plot viewer", () => addViewer(page, "3d scatter plot"));
    await session.step(16, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
    await run.scenario("The axes are assigned automatically", async () => {
      await session.step(19, "Then X column input in 3d scatter plot viewer should contain text \"AGE\"", () => shouldContainText(page, el("X column input in 3d scatter plot viewer"), "AGE"));
      await session.step(20, "And Y column input in 3d scatter plot viewer should contain text \"HEIGHT\"", () => shouldContainText(page, el("Y column input in 3d scatter plot viewer"), "HEIGHT"));
      await session.step(21, "And Z column input in 3d scatter plot viewer should contain text \"WEIGHT\"", () => shouldContainText(page, el("Z column input in 3d scatter plot viewer"), "WEIGHT"));
      await session.step(22, "And properties of 3d scatter plot viewer should be:", () => propertiesShouldBe(page, el("3d scatter plot viewer"), [["X","AGE"],["Y","HEIGHT"],["Z","WEIGHT"]]));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reassigning X and Z moves the selectors and redraws the scene", async () => {
      await session.step(29, "When user sets properties of 3d scatter plot viewer:", () => setProperties(page, el("3d scatter plot viewer"), [["X","WEIGHT"],["Z","AGE"]]));
      await session.step(32, "Then X column input in 3d scatter plot viewer should contain text \"WEIGHT\"", () => shouldContainText(page, el("X column input in 3d scatter plot viewer"), "WEIGHT"));
      await session.step(33, "And Z column input in 3d scatter plot viewer should contain text \"AGE\"", () => shouldContainText(page, el("Z column input in 3d scatter plot viewer"), "AGE"));
      await session.step(34, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(35, "When user sets properties of 3d scatter plot viewer:", () => setProperties(page, el("3d scatter plot viewer"), [["X","AGE"],["Z","WEIGHT"]]));
      await session.step(38, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color by a category draws a legend of its values", async () => {
      await session.step(42, "Then legend of 3d scatter plot viewer should be hidden", () => shouldBe(page, el("legend of 3d scatter plot viewer"), "hidden"));
      await session.step(43, "When user sets \"Color\" property of 3d scatter plot viewer to \"SEX\"", () => setProperty(page, "Color", el("3d scatter plot viewer"), "SEX"));
      await session.step(44, "Then legend of 3d scatter plot viewer should have 2 items", () => shouldHaveItems(page, el("legend of 3d scatter plot viewer"), 2));
      await session.step(45, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color by a number takes the legend away", async () => {
      await session.step(49, "When user sets \"Color\" property of 3d scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("3d scatter plot viewer"), "AGE"));
      await session.step(50, "Then legend of 3d scatter plot viewer should be hidden", () => shouldBe(page, el("legend of 3d scatter plot viewer"), "hidden"));
      await session.step(51, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Marker type redraws the markers", async () => {
      await session.step(55, "When user sets \"Marker Type\" property of 3d scatter plot viewer to \"box\"", () => setProperty(page, "Marker Type", el("3d scatter plot viewer"), "box"));
      await session.step(56, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(57, "When user sets \"Marker Type\" property of 3d scatter plot viewer to \"sphere\"", () => setProperty(page, "Marker Type", el("3d scatter plot viewer"), "sphere"));
      await session.step(58, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(59, "When user sets \"Marker Type\" property of 3d scatter plot viewer to \"cylinder\"", () => setProperty(page, "Marker Type", el("3d scatter plot viewer"), "cylinder"));
      await session.step(60, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(61, "When user sets \"Marker Type\" property of 3d scatter plot viewer to \"octahedron\"", () => setProperty(page, "Marker Type", el("3d scatter plot viewer"), "octahedron"));
      await session.step(62, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Marker opacity redraws the markers", async () => {
      await session.step(66, "When user sets \"Marker Opacity\" property of 3d scatter plot viewer to \"25\"", () => setProperty(page, "Marker Opacity", el("3d scatter plot viewer"), "25"));
      await session.step(67, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(68, "When user sets \"Marker Opacity\" property of 3d scatter plot viewer to \"69\"", () => setProperty(page, "Marker Opacity", el("3d scatter plot viewer"), "69"));
      await session.step(69, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Axes hides and restores the axes", async () => {
      await session.step(73, "When user sets \"Show Axes\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Axes", el("3d scatter plot viewer"), "false"));
      await session.step(74, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(75, "When user sets \"Show Axes\" property of 3d scatter plot viewer to \"true\"", () => setProperty(page, "Show Axes", el("3d scatter plot viewer"), "true"));
      await session.step(76, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic X axis re-scales the scene without error", async () => {
      await session.step(80, "When user sets \"X Axis Type\" property of 3d scatter plot viewer to \"logarithmic\"", () => setProperty(page, "X Axis Type", el("3d scatter plot viewer"), "logarithmic"));
      await session.step(81, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
      await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(84, "When user sets \"X Axis Type\" property of 3d scatter plot viewer to \"linear\"", () => setProperty(page, "X Axis Type", el("3d scatter plot viewer"), "linear"));
      await session.step(85, "Then the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A drag rotates the scene and Reset View brings the camera home", async () => {
      await session.step(89, "When user drags across the \"view\" area of 3d scatter plot viewer", () => dragAcrossArea(page, "view", el("3d scatter plot viewer")));
      await session.step(90, "Then the \"camera x\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "camera x", el("3d scatter plot viewer")));
      await session.step(91, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(92, "When user picks \"Reset View\" from the context menu of 3d scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("3d scatter plot viewer")));
      await session.step(93, "Then the \"camera x\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "camera x", el("3d scatter plot viewer"), 0));
      await session.step(94, "And the \"camera y\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "camera y", el("3d scatter plot viewer"), 0));
      await session.step(95, "And the \"camera distance\" reading of 3d scatter plot viewer should be 4", () => readingIs(page, "camera distance", el("3d scatter plot viewer"), 4));
      await session.step(96, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The mouse wheel zooms the scene", async () => {
      await session.step(100, "When user scrolls the mouse wheel up over the \"view\" area of 3d scatter plot viewer", () => wheelOverArea(page, "up", "view", el("3d scatter plot viewer")));
      await session.step(101, "Then the \"camera distance\" reading of 3d scatter plot viewer should be lower than before", () => readingLower(page, "camera distance", el("3d scatter plot viewer")));
      await session.step(102, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(103, "When user scrolls the mouse wheel down over the \"view\" area of 3d scatter plot viewer", () => wheelOverArea(page, "down", "view", el("3d scatter plot viewer")));
      await session.step(104, "Then the \"camera distance\" reading of 3d scatter plot viewer should be higher than before", () => readingHigher(page, "camera distance", el("3d scatter plot viewer")));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click makes a row current, Shift-click selects it", async () => {
      await session.step(108, "When user clicks on the \"point\" area of 3d scatter plot viewer", () => clickArea(page, "point", el("3d scatter plot viewer")));
      await session.step(109, "Then the \"current row\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "current row", el("3d scatter plot viewer")));
      await session.step(110, "And the table should have a current row", () => hasCurrentRow(page));
      await session.step(111, "When user clears the row selection", () => clearSelection(page));
      await session.step(112, "And user clicks on the \"point\" area of 3d scatter plot viewer holding Shift", () => clickAreaHolding(page, "point", el("3d scatter plot viewer"), "Shift"));
      await session.step(113, "Then some rows should be selected", () => someSelected(page));
      await session.step(114, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(115, "When user clears the row selection", () => clearSelection(page));
      await session.step(116, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Filtered Out Points brings the filtered-away rows back", async () => {
      await session.step(119, "When user filters rows where \"SEX\" is \"F\"", () => filterTo(page, "SEX", "F"));
      await session.step(120, "Then 553 rows should pass the filter", () => filterPasses(page, 553));
      await session.step(121, "And the \"rows shown\" reading of 3d scatter plot viewer should be lower than before", () => readingLower(page, "rows shown", el("3d scatter plot viewer")));
      await session.step(122, "When user sets \"Show Filtered Out Points\" property of 3d scatter plot viewer to \"true\"", () => setProperty(page, "Show Filtered Out Points", el("3d scatter plot viewer"), "true"));
      await session.step(123, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(124, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(125, "When user sets \"Show Filtered Out Points\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Filtered Out Points", el("3d scatter plot viewer"), "false"));
      await session.step(126, "Then the \"rows shown\" reading of 3d scatter plot viewer should be lower than before", () => readingLower(page, "rows shown", el("3d scatter plot viewer")));
      await session.step(127, "When user resets the filter", () => resetFilter(page));
      await session.step(128, "Then the \"rows shown\" reading of 3d scatter plot viewer should be 872", () => readingIs(page, "rows shown", el("3d scatter plot viewer"), 872));
      await session.step(129, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A hover on a bar chart highlights the matching points", async () => {
      await session.step(132, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","SEX"]]));
      await session.step(134, "Then bar chart viewer should have a \"bar F\" area", () => hasArea(page, el("bar chart viewer"), "bar F"));
      await session.step(135, "And \"Show Mouse Over Row Group\" property of 3d scatter plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Mouse Over Row Group", el("3d scatter plot viewer"), "true"));
      await session.step(136, "And the \"highlighted rows\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "highlighted rows", el("3d scatter plot viewer"), 0));
      await session.step(137, "When user takes a snapshot of 3d scatter plot viewer", () => takeSnapshot(page, el("3d scatter plot viewer")));
      await session.step(138, "And user hovers over the \"bar F\" area of bar chart viewer", () => hoverArea(page, "bar F", el("bar chart viewer")));
      await session.step(139, "Then the \"highlighted rows\" reading of 3d scatter plot viewer should be higher than before", () => readingHigher(page, "highlighted rows", el("3d scatter plot viewer")));
      await session.step(140, "And the \"scene signature\" reading of 3d scatter plot viewer should differ from before", () => readingDiffers(page, "scene signature", el("3d scatter plot viewer")));
      await session.step(141, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(142, "And user sets \"Show Mouse Over Row Group\" property of 3d scatter plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("3d scatter plot viewer"), "false"));
      await session.step(143, "And user hovers over the \"bar F\" area of bar chart viewer", () => hoverArea(page, "bar F", el("bar chart viewer")));
      await session.step(144, "Then the \"highlighted rows\" reading of 3d scatter plot viewer should be 0", () => readingIs(page, "highlighted rows", el("3d scatter plot viewer"), 0));
      await session.step(145, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
      await session.step(146, "And user sets \"Show Mouse Over Row Group\" property of 3d scatter plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("3d scatter plot viewer"), "true"));
      await session.step(147, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(148, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position docks the legend on the side it names", async () => {
      await session.step(152, "When user sets properties of 3d scatter plot viewer:", () => setProperties(page, el("3d scatter plot viewer"), [["Color","SEX"],["Legend Visibility","Always"]]));
      await session.step(155, "Then legend of 3d scatter plot viewer should have 2 items", () => shouldHaveItems(page, el("legend of 3d scatter plot viewer"), 2));
      await session.step(156, "When user sets \"Legend Position\" property of 3d scatter plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("3d scatter plot viewer"), "Left"));
      await session.step(157, "Then the legend of 3d scatter plot viewer should be on the left", () => legendSide(page, el("3d scatter plot viewer"), "left"));
      await session.step(158, "When user sets \"Legend Position\" property of 3d scatter plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("3d scatter plot viewer"), "Right"));
      await session.step(159, "Then the legend of 3d scatter plot viewer should be on the right", () => legendSide(page, el("3d scatter plot viewer"), "right"));
      await session.step(160, "When user sets properties of 3d scatter plot viewer:", () => setProperties(page, el("3d scatter plot viewer"), [["Legend Position","Auto"],["Legend Visibility","Auto"],["Color",""]]));
      await session.step(164, "Then legend of 3d scatter plot viewer should be hidden", () => shouldBe(page, el("legend of 3d scatter plot viewer"), "hidden"));
      await session.step(165, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
