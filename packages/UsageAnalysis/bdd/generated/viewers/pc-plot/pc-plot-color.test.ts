/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-color.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {colorAgain, colorCodedAs, colorCodedCategorically, colorConditional, colorOff, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaPainted, areaRepainted, hasArea, hasNoArea, legendItemsDiffer, legendLists, legendSameItems, legendSide, noErrors, painted, pickFromContextMenu, propertyShouldBe, repaintedBy, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot colouring, legend and colour scale", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-color.feature", import.meta.url);
  test("PC plot colouring, legend and colour scale", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(19, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(20, "And \"Color\" property of pc plot viewer should be \"\"", () => propertyShouldBe(page, "Color", el("pc plot viewer"), ""));
    await session.step(21, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
    await session.step(22, "And pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
    await run.scenario("A categorical colour column lists its categories in the legend", async () => {
      await session.step(25, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(26, "Then \"Color\" property of pc plot viewer should be \"RACE\"", () => propertyShouldBe(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(27, "And legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(28, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(29, "And \"Caucasian\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Caucasian\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(30, "And \"Asian\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Asian\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(31, "And \"Black\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Black\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(32, "And \"Other\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Other\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(33, "And the \"Caucasian\" and \"Asian\" items in the legend of pc plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("pc plot viewer")));
      await session.step(34, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(35, "And pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(36, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(37, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(38, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Position moves the legend and Legend Visibility takes it away", async () => {
      await session.step(42, "When user sets \"Color\" property of pc plot viewer to \"RACE\"", () => setProperty(page, "Color", el("pc plot viewer"), "RACE"));
      await session.step(43, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(44, "When user sets \"Legend Position\" property of pc plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Left"));
      await session.step(45, "Then the legend of pc plot viewer should be on the left", () => legendSide(page, el("pc plot viewer"), "left"));
      await session.step(46, "And the legend of pc plot viewer should list the same items as before", () => legendSameItems(page, el("pc plot viewer")));
      await session.step(47, "When user sets \"Legend Position\" property of pc plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Right"));
      await session.step(48, "Then the legend of pc plot viewer should be on the right", () => legendSide(page, el("pc plot viewer"), "right"));
      await session.step(49, "And the legend of pc plot viewer should list the same items as before", () => legendSameItems(page, el("pc plot viewer")));
      await session.step(50, "When user sets \"Legend Position\" property of pc plot viewer to \"Top\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Top"));
      await session.step(51, "Then the legend of pc plot viewer should be on the top", () => legendSide(page, el("pc plot viewer"), "top"));
      await session.step(52, "When user sets \"Legend Position\" property of pc plot viewer to \"Bottom\"", () => setProperty(page, "Legend Position", el("pc plot viewer"), "Bottom"));
      await session.step(53, "Then \"Legend Position\" property of pc plot viewer should be \"Bottom\"", () => propertyShouldBe(page, "Legend Position", el("pc plot viewer"), "Bottom"));
      await session.step(54, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(55, "When user sets \"Legend Visibility\" property of pc plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("pc plot viewer"), "Never"));
      await session.step(56, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(57, "When user sets \"Legend Visibility\" property of pc plot viewer to \"Auto\"", () => setProperty(page, "Legend Visibility", el("pc plot viewer"), "Auto"));
      await session.step(58, "Then legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(59, "And the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(60, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Legend Position","Auto"],["Color",""]]), [["Legend Position","Auto"],["Color",""]]);
      await session.step(63, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numerical colour column draws the colour scale instead of a legend", async () => {
      await session.step(67, "When user sets \"Color\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Color", el("pc plot viewer"), "AGE"));
      await session.step(68, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(69, "And pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(70, "And the \"color scale\" area of pc plot viewer should be painted", () => areaPainted(page, "color scale", el("pc plot viewer")));
      await session.step(71, "And the \"color scale\" area of pc plot viewer should contain the color \"#FF0000\"", () => areaColor(page, "color scale", el("pc plot viewer"), "#FF0000"));
      await session.step(72, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(73, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(74, "Then pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Scheme > Invert Color Scheme repaints the scale and the lines", async () => {
      await session.step(78, "When user sets \"Color\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Color", el("pc plot viewer"), "AGE"));
      await session.step(79, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "false"));
      await session.step(80, "When user picks \"Color Scheme > Invert Color Scheme\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Invert Color Scheme", el("pc plot viewer")));
      await session.step(81, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "true"));
      await session.step(82, "And the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(83, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(84, "When user picks \"Color Scheme > Invert Color Scheme\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Invert Color Scheme", el("pc plot viewer")));
      await session.step(85, "Then \"Invert Color Scheme\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Invert Color Scheme", el("pc plot viewer"), "false"));
      await session.step(86, "And the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(87, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Scheme > Edit... opens the column's colour-coding dialog", async () => {
      await session.step(91, "When user picks \"Color Scheme > Edit...\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Color Scheme > Edit...", el("pc plot viewer")));
      await session.step(92, "Then \"Color-coding: AGE\" dialog should be visible", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "visible"));
      await session.step(93, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(94, "Then \"Color-coding: AGE\" dialog should be absent", () => shouldBe(page, el("\"Color-coding: AGE\" dialog"), "absent"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Color Min and Color Max clamp the scale and recolour the lines", async () => {
      await session.step(98, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Min","30"],["Color Max","60"]]), [["Color Min","30"],["Color Max","60"]]);
      await session.step(101, "Then the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(102, "And pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(103, "When user sets \"Color Axis Type\" property of pc plot viewer to \"logarithmic\"", () => setProperty(page, "Color Axis Type", el("pc plot viewer"), "logarithmic"));
      await session.step(104, "Then the \"color scale\" area of pc plot viewer should have repainted", () => areaRepainted(page, "color scale", el("pc plot viewer")));
      await session.step(105, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Axis Type","linear"],["Color Min",""],["Color Max",""]]), [["Color Axis Type","linear"],["Color Min",""],["Color Max",""]]);
      await session.step(109, "Then pc plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("pc plot viewer"), 2000));
      await session.step(110, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(111, "Then pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A conditional colour coding on the column hands the plot its bins", async () => {
      await session.step(115, "When user sets \"Color\" property of pc plot viewer to \"HEIGHT\"", () => setProperty(page, "Color", el("pc plot viewer"), "HEIGHT"));
      await session.step(116, "Then pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(117, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["20-150","#00FF00"],["150-250","#FFA500"]]), [["20-150","#00FF00"],["150-250","#FFA500"]]);
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
    await run.scenario("HEIGHT's colour coding switched categorical, numerical and off leaves no legend behind", async () => {
      await session.step(134, "When user sets \"Color\" property of pc plot viewer to \"HEIGHT\"", () => setProperty(page, "Color", el("pc plot viewer"), "HEIGHT"));
      await session.step(135, "And user colors \"HEIGHT\" column categorically again", () => colorAgain(page, "HEIGHT", "categorically"));
      await session.step(136, "Then \"HEIGHT\" column should be color-coded categorically", () => colorCodedCategorically(page, "HEIGHT"));
      await session.step(137, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(138, "When user colors \"HEIGHT\" column linearly again", () => colorAgain(page, "HEIGHT", "linearly"));
      await session.step(139, "Then \"HEIGHT\" column should be color-coded linearly", () => colorCodedAs(page, "HEIGHT", "linearly"));
      await session.step(140, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(141, "And pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(142, "When user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
      await session.step(143, "Then \"HEIGHT\" column should have no color coding", () => noColorCoding(page, "HEIGHT"));
      await session.step(144, "And legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(145, "And pc plot viewer should have a \"color scale\" area", () => hasArea(page, el("pc plot viewer"), "color scale"));
      await session.step(146, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(147, "When user sets \"Color\" property of pc plot viewer to \"\"", () => setProperty(page, "Color", el("pc plot viewer"), ""));
      await session.step(148, "Then pc plot viewer should not have a \"color scale\" area", () => hasNoArea(page, el("pc plot viewer"), "color scale"));
      await session.step(149, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A DateTime colour column is split by its Color Map", async () => {
      await session.step(152, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color","STARTED"],["Color Map","year"]]), [["Color","STARTED"],["Color Map","year"]]);
      await session.step(155, "Then legend of pc plot viewer should be visible", () => shouldBe(page, el("legend of pc plot viewer"), "visible"));
      await session.step(156, "And the legend of pc plot viewer should list 3 items", () => legendLists(page, el("pc plot viewer"), 3));
      await session.step(157, "And \"1990\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"1990\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(158, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(159, "When user sets \"Color Map\" property of pc plot viewer to \"quarter\"", () => setProperty(page, "Color Map", el("pc plot viewer"), "quarter"));
      await session.step(160, "Then the legend of pc plot viewer should list 4 items", () => legendLists(page, el("pc plot viewer"), 4));
      await session.step(161, "And \"Q1\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"Q1\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(162, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(163, "When user sets \"Color Map\" property of pc plot viewer to \"month\"", () => setProperty(page, "Color Map", el("pc plot viewer"), "month"));
      await session.step(164, "Then the legend of pc plot viewer should list 12 items", () => legendLists(page, el("pc plot viewer"), 12));
      await session.step(165, "And \"January\" legend item in legend of pc plot viewer should be visible", () => shouldBe(page, el("\"January\" legend item in legend of pc plot viewer"), "visible"));
      await session.step(166, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Color Map",""],["Color",""]]), [["Color Map",""],["Color",""]]);
      await session.step(169, "Then legend of pc plot viewer should be hidden", () => shouldBe(page, el("legend of pc plot viewer"), "hidden"));
      await session.step(170, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
      await session.step(171, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
