/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-color.feature
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorCodedAs, colorConditional, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaPainted, areaRepainted, hasArea, hasNoArea, legendItemsDiffer, legendLists, legendSameItems, legendSide, noErrors, pickFromContextMenu, propertyShouldBe, repaintedBy, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot colouring, legend and colour scale", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-color.feature", import.meta.url);
  test("PC plot colouring, legend and colour scale", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
    await session.step(18, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And \"Color\" property of pc plot viewer should be \"\"", () => propertyShouldBe(page, "Color", el("pc plot viewer"), ""));
    await session.step(20, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
    await session.step(21, "And pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
    await run.scenario("A categorical colour column lists its categories in the legend", async () => {
      await session.step(24, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(25, "Then \"Color\" property of pc plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(26, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(27, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(28, "And \"Caucasian\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(29, "And \"Asian\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(30, "And \"Black\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(31, "And \"Other\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Other\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(32, "And the \"Caucasian\" and \"Asian\" items in the legend of pc plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("pc plot viewer")));
      await session.step(33, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(34, "And pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(35, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(36, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(37, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position moves the legend and Legend Visibility takes it away", async () => {
      await session.step(41, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(42, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(43, "When user sets \"Legend Position\" property of pc plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Left"));
      await session.step(44, "Then the legend of pc plot viewer should be on the left", () => legendSide(page, el("pc plot viewer"), "left"));
      await session.step(45, "And the legend of pc plot viewer should list the same items as before", () => legendSameItems(page, el("pc plot viewer")));
      await session.step(46, "When user sets \"Legend Position\" property of pc plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Right"));
      await session.step(47, "Then the legend of pc plot viewer should be on the right", () => legendSide(page, el("pc plot viewer"), "right"));
      await session.step(48, "And the legend of pc plot viewer should list the same items as before", () => legendSameItems(page, el("pc plot viewer")));
      await session.step(49, "When user sets \"Legend Position\" property of pc plot viewer to \"Top\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Top"));
      await session.step(50, "Then the legend of pc plot viewer should be on the top", () => legendSide(page, el("pc plot viewer"), "top"));
      await session.step(51, "When user sets \"Legend Position\" property of pc plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Bottom"));
      await session.step(52, "Then \"Legend Position\" property of pc plot viewer should be \"Bottom\"", () => propertyShouldBe(page, "Legend Position", el("pc plot viewer"), "Bottom"));
      await session.step(53, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(54, "When user sets \"Legend Visibility\" property of pc plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("pc plot viewer"), "Never"));
      await session.step(55, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(56, "When user sets \"Legend Visibility\" property of pc plot viewer to \"Auto\"", () => setProperty(page, "Legend Visibility", el("pc plot viewer"), "Auto"));
      await session.step(57, "Then legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(58, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(59, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Legend Position","Auto"],["Color",""]]));
      await session.step(62, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical colour column draws the colour scale instead of a legend", async () => {
      await session.step(66, "When user sets \"Color\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Color", el("pc plot viewer"), "AGE"));
      await session.step(67, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(68, "And pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(69, "And the \"color scale\" area of pc plot viewer should be painted", () => areaPainted(page, "color scale", el("pc plot viewer")));
      await session.step(70, "And the \"color scale\" area of pc plot viewer should contain the color \"#FF0000\"", () => areaColor(page, "color scale", el("pc plot viewer"), "#FF0000"));
      await session.step(71, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(72, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(73, "Then pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Scheme > Invert Color Scheme repaints the scale and the lines", async () => {
      await session.step(77, "When user sets \"Color\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Color", el("pc plot viewer"), "AGE"));
      await session.step(78, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "false"));
      await session.step(79, "When user picks \"Color Scheme > Invert Color Scheme\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Invert Color Scheme", el("pc plot viewer")));
      await session.step(80, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "true"));
      await session.step(81, "And the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(82, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(83, "When user picks \"Color Scheme > Invert Color Scheme\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Invert Color Scheme", el("pc plot viewer")));
      await session.step(84, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "false"));
      await session.step(85, "And the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(86, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Scheme > Edit... opens the column's colour-coding dialog", async () => {
      await session.step(90, "When user picks \"Color Scheme > Edit...\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Edit...", el("pc plot viewer")));
      await session.step(91, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(92, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(93, "Then \"Color-coding: AGE\" dialog should be absent", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "absent"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Min and Color Max clamp the scale and recolour the lines", async () => {
      await session.step(97, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Min","30"],["Color Max","60"]]));
      await session.step(100, "Then \"Color Min\" property of pc plot viewer should be \"30\"", () => propertyShouldBe(page, "Color Min", el("pc plot viewer"), "30"));
      await session.step(101, "And the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(102, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(103, "When user sets \"Color Axis Type\" property of pc plot viewer to \"logarithmic\"", () => setProperty(page, "Color Axis Type", el("pc plot viewer"), "logarithmic"));
      await session.step(104, "Then the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(105, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Axis Type","linear"],["Color Min",""],["Color Max",""]]));
      await session.step(109, "Then pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(110, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(111, "Then pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A conditional colour coding on the column hands the plot its bins", async () => {
      await session.step(115, "When user sets \"Color\" property of pc plot viewer to \"HEIGHT\"", () => setProperty(page, "Color", el("pc plot viewer"), "HEIGHT"));
      await session.step(116, "Then pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(117, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["20-150","#00FF00"],["150-250","#FFA500"]]));
      await session.step(120, "Then \"HEIGHT\" column should be color-coded conditionally", () => colorCodedAs(page, "HEIGHT", "conditionally"));
      await session.step(121, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(122, "And the legend of pc plot viewer should list 3 items", () => legendLists(page, el("pc plot viewer"), 3));
      await session.step(123, "And \"20-150\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"20-150\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(124, "And \"150-250\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"150-250\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(125, "And pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(126, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(127, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(128, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(129, "And pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(130, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(131, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A DateTime colour column is split by its Color Map", async () => {
      await session.step(134, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color","STARTED"],["Color Map","year"]]));
      await session.step(137, "Then legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(138, "And the legend of pc plot viewer should list 3 items", () => legendLists(page, el("pc plot viewer"), 3));
      await session.step(139, "And \"1990\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"1990\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(140, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(141, "When user sets \"Color Map\" property of pc plot viewer to \"quarter\"", () => setProperty(page, "Color Map", el("pc plot viewer"), "quarter"));
      await session.step(142, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(143, "And \"Q1\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Q1\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(144, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(145, "When user sets \"Color Map\" property of pc plot viewer to \"month\"", () => setProperty(page, "Color Map", el("pc plot viewer"), "month"));
      await session.step(146, "Then the legend of pc plot viewer should list 12 items", () => legendLists(page, el("pc plot viewer"), 12));
      await session.step(147, "And \"January\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"January\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(148, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Map",""],["Color",""]]));
      await session.step(151, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(152, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(153, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
