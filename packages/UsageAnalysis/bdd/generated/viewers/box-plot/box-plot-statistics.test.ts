/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-statistics.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, colorCategorical, colorConditional, colorLinear, colorLinearOver, colorOff, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaLessInk, areaMoreInk, areaRepainted, areasDiffer, closeContextMenu, dragSelectionOverArea, hasArea, hasNoArea, hoverArea, moreHighlight, noBalloons, noErrors, pickFromAreaContextMenu, pointerAway, propertiesShouldBe, propertyShouldBe, repainted, setProperties, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot statistics and coloring", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-statistics.feature", import.meta.url);
  test("Box plot statistics and coloring", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]), [["Value","AGE"],["Category 1","SEX"]]);
    await run.scenario("Box coloring", async () => {
      await session.step(17, "Then \"Whisker Color\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Whisker Color", el("box plot viewer"), ""));
      await session.step(18, "And the \"M values\" and \"F values\" areas of box plot viewer should be painted in different colors", () => areasDiffer(page, "M values", "F values", el("box plot viewer")));
      await session.step(19, "When user sets \"Whisker Color\" property of box plot viewer to \"#1F77B4\"", () => setProperty(page, "Whisker Color", el("box plot viewer"), "#1F77B4"));
      await session.step(20, "Then \"Whisker Color\" property of box plot viewer should be \"#1F77B4\"", () => propertyShouldBe(page, "Whisker Color", el("box plot viewer"), "#1F77B4"));
      await session.step(21, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(22, "And the \"M values\" area of box plot viewer should contain the color \"#1F77B4\"", () => areaColor(page, "M values", el("box plot viewer"), "#1F77B4"));
      await session.step(23, "And the \"F values\" area of box plot viewer should contain the color \"#1F77B4\"", () => areaColor(page, "F values", el("box plot viewer"), "#1F77B4"));
      await session.step(24, "When user sets \"Whisker Color\" property of box plot viewer to \"\"", () => setProperty(page, "Whisker Color", el("box plot viewer"), ""));
      await session.step(25, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
    });
    await run.scenario("The statistics strip and its ladder", async () => {
      await session.step(28, "Then \"Show Statistics\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Statistics", el("box plot viewer"), "true"));
      await session.step(29, "And box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
      await session.step(30, "When user hovers over the \"stats\" area of box plot viewer", () => hoverArea(page, "stats", el("box plot viewer")));
      await session.step(31, "Then close-stats icon in box plot viewer should be visible", () => shouldBe(page, el("close-stats icon in box plot viewer"), "visible"));
      await session.step(32, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(33, "Then close-stats icon in box plot viewer should be hidden", () => shouldBe(page, el("close-stats icon in box plot viewer"), "hidden"));
      await session.step(34, "And \"Show Avg\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Avg", el("box plot viewer"), "true"));
      await session.step(35, "And \"Show Med\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Med", el("box plot viewer"), "true"));
      await session.step(36, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]);
      await session.step(43, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(44, "And box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
      await session.step(45, "And the \"stats\" area of box plot viewer should have more ink than before", () => areaMoreInk(page, "stats", el("box plot viewer")));
      await session.step(46, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]);
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
      await session.step(54, "When user sets \"Statistics Format\" property of box plot viewer to \"#,##0.00\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "#,##0.00"));
      await session.step(55, "Then the \"stats\" area of box plot viewer should have repainted", () => areaRepainted(page, "stats", el("box plot viewer")));
      await session.step(56, "And \"Statistics Format\" property of box plot viewer should be \"#,##0.00\"", () => propertyShouldBe(page, "Statistics Format", el("box plot viewer"), "#,##0.00"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
      await session.step(58, "When user sets \"Statistics Format\" property of box plot viewer to \"auto\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "auto"));
      await session.step(59, "And user picks \"Show Total Count\" from the context menu of the \"stats\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Show Total Count", "stats", el("box plot viewer")));
      await session.step(60, "Then \"Show Total Count\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Total Count", el("box plot viewer"), "false"));
      await session.step(61, "When user picks \"Show Total Count\" from the context menu of the \"stats\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Show Total Count", "stats", el("box plot viewer")));
      await session.step(62, "Then \"Show Total Count\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Total Count", el("box plot viewer"), "true"));
      await session.step(63, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(64, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Total Count","false"],["Show Inliers Count","false"],["Show Outliers Count","false"],["Show Stdev","false"],["Show Q1","false"],["Show Q3","false"]]), [["Show Total Count","false"],["Show Inliers Count","false"],["Show Outliers Count","false"],["Show Stdev","false"],["Show Q1","false"],["Show Q3","false"]]);
    });
    await run.scenario("The T key toggles the p-value", async () => {
      await session.step(73, "Then box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
      await session.step(74, "When user sets \"Show P Value\" property of box plot viewer to \"false\"", () => setProperty(page, "Show P Value", el("box plot viewer"), "false"));
      await session.step(75, "Then box plot viewer should not have a \"p value\" area", () => hasNoArea(page, el("box plot viewer"), "p value"));
      await session.step(76, "When user presses t in box plot viewer", () => pressKeyIn(page, "t", el("box plot viewer")));
      await session.step(77, "Then \"Show P Value\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "true"));
      await session.step(78, "And box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
      await session.step(79, "When user presses t in box plot viewer", () => pressKeyIn(page, "t", el("box plot viewer")));
      await session.step(80, "Then \"Show P Value\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "false"));
      await session.step(81, "And box plot viewer should not have a \"p value\" area", () => hasNoArea(page, el("box plot viewer"), "p value"));
      await session.step(82, "When user sets \"Show P Value\" property of box plot viewer to \"true\"", () => setProperty(page, "Show P Value", el("box plot viewer"), "true"));
    });
    await run.scenario("Three groups take the Alexander-Govern branch", async () => {
      await session.step(85, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(86, "Then box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
      await session.step(87, "And box plot viewer should have a \"category Asian\" area", () => hasArea(page, el("box plot viewer"), "category Asian"));
      await session.step(88, "And box plot viewer should have a \"category Black\" area", () => hasArea(page, el("box plot viewer"), "category Black"));
      await session.step(89, "And box plot viewer should have a \"category Caucasian\" area", () => hasArea(page, el("box plot viewer"), "category Caucasian"));
      await session.step(90, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(91, "Then tooltip should contain text \"Alexander\"", () => shouldContainText(page, el("tooltip"), "Alexander"));
      await session.step(92, "And show-group-stats icon in box plot viewer should be visible", () => shouldBe(page, el("show-group-stats icon in box plot viewer"), "visible"));
      await session.step(93, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(94, "Then show-group-stats icon in box plot viewer should be hidden", () => shouldBe(page, el("show-group-stats icon in box plot viewer"), "hidden"));
      await session.step(95, "When user sets \"Category 1\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 1", el("box plot viewer"), "SEX"));
    });
    await run.scenario("The violin style", async () => {
      await session.step(98, "When user sets \"Plot Style\" property of box plot viewer to \"violin\"", () => setProperty(page, "Plot Style", el("box plot viewer"), "violin"));
      await session.step(99, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(100, "And the \"M values\" area of box plot viewer should have more ink than before", () => areaMoreInk(page, "M values", el("box plot viewer")));
      await session.step(101, "And the \"F values\" area of box plot viewer should have more ink than before", () => areaMoreInk(page, "F values", el("box plot viewer")));
      await session.step(102, "When user sets \"Bins\" property of box plot viewer to \"50\"", () => setProperty(page, "Bins", el("box plot viewer"), "50"));
      await session.step(103, "And user sets \"Bins\" property of box plot viewer to \"500\"", () => setProperty(page, "Bins", el("box plot viewer"), "500"));
      await session.step(104, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(105, "When user sets \"Interquartile Line Width\" property of box plot viewer to \"10\"", () => setProperty(page, "Interquartile Line Width", el("box plot viewer"), "10"));
      await session.step(106, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(107, "When user sets \"Violin Line Width\" property of box plot viewer to \"4\"", () => setProperty(page, "Violin Line Width", el("box plot viewer"), "4"));
      await session.step(108, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(109, "When user sets \"Violin Whisker Color\" property of box plot viewer to \"#00AA00\"", () => setProperty(page, "Violin Whisker Color", el("box plot viewer"), "#00AA00"));
      await session.step(110, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(111, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Plot Style","box"],["Bins","100"],["Interquartile Line Width","6"],["Violin Line Width","2"]]), [["Plot Style","box"],["Bins","100"],["Interquartile Line Width","6"],["Violin Line Width","2"]]);
      await session.step(116, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(117, "And the \"M values\" area of box plot viewer should have less ink than before", () => areaLessInk(page, "M values", el("box plot viewer")));
    });
    await run.scenario("Column color coding drives the marker colors", async () => {
      await session.step(120, "When user colors \"WEIGHT\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "WEIGHT", "#0000FF", "#FF0000"));
      await session.step(121, "And user sets \"Marker Color Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "WEIGHT"));
      await session.step(122, "Then \"Marker Color Column\" property of box plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "WEIGHT"));
      await session.step(123, "And box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
      await session.step(124, "When user colors \"WEIGHT\" column linearly from \"#0000FF\" to \"#FF0000\" over 60 to 120", () => colorLinearOver(page, "WEIGHT", "#0000FF", "#FF0000", 60, 120));
      await session.step(125, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(126, "When user clears the row selection", () => clearSelection(page));
      await session.step(127, "And user drags a selection box over the \"M values\" area of box plot viewer", () => dragSelectionOverArea(page, "M values", el("box plot viewer")));
      await session.step(128, "Then some rows should be selected", () => someSelected(page));
      await session.step(129, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(130, "When user colors \"WEIGHT\" column conditionally:", () => colorConditional(page, "WEIGHT", [["50-90","#00FF00"],["90-150","#800080"]]), [["50-90","#00FF00"],["90-150","#800080"]]);
      await session.step(133, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(134, "And the \"F values\" area of box plot viewer should contain the color \"#00FF00\"", () => areaColor(page, "F values", el("box plot viewer"), "#00FF00"));
      await session.step(135, "And box plot viewer should show a selection highlight", () => someHighlight(page, el("box plot viewer")));
      await session.step(136, "When user clears the row selection", () => clearSelection(page));
      await session.step(137, "And user removes the coloring of \"WEIGHT\" column", () => colorOff(page, "WEIGHT"));
      await session.step(138, "And user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(139, "And user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#E41A1C"]]), [["M","#E41A1C"]]);
      await session.step(141, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(142, "And the \"M values\" area of box plot viewer should contain the color \"#E41A1C\"", () => areaColor(page, "M values", el("box plot viewer"), "#E41A1C"));
      await session.step(143, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(144, "And user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("A datetime value", async () => {
      await session.step(147, "When user sets \"Value\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Value", el("box plot viewer"), "STARTED"));
      await session.step(148, "Then \"Value\" property of box plot viewer should be \"STARTED\"", () => propertyShouldBe(page, "Value", el("box plot viewer"), "STARTED"));
      await session.step(149, "And box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
      await session.step(150, "And no errors should have been logged", () => noErrors(page));
      await session.step(151, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(152, "When user sets \"Statistics Format\" property of box plot viewer to \"yyyy-MM-dd HH:mm\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "yyyy-MM-dd HH:mm"));
      await session.step(153, "Then \"Statistics Format\" property of box plot viewer should be \"yyyy-MM-dd HH:mm\"", () => propertyShouldBe(page, "Statistics Format", el("box plot viewer"), "yyyy-MM-dd HH:mm"));
      await session.step(154, "And box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
      await session.step(155, "And no errors should have been logged", () => noErrors(page));
      await session.step(156, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(157, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","AGE"],["Statistics Format","auto"]]), [["Value","AGE"],["Statistics Format","auto"]]);
    });
    run.finish();
  });
});
