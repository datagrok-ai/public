/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot-statistics.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clickEmptySpace} from '../../bindings/box-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, colorCategorical, colorConditional, colorLinear, colorLinearOver, colorOff, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, dragSelectionOverArea, hasArea, hoverArea, moreHighlight, noErrors, painted, paintedInColors, pickFromAreaContextMenu, pointerAway, propertiesShouldBe, propertyShouldBe, repainted, setProperties, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot statistics and coloring", () => {
  const session = feature(test, "features/viewers/box-plot-statistics.feature", import.meta.url);
  test("Box plot statistics and coloring", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]));
    await run.scenario("Box coloring", async () => {
      await session.step(16, "Then \"Whisker Color\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Whisker Color", el("box plot viewer"), ""));
      await session.step(17, "And box plot viewer should be painted in at least 2 colors", () => paintedInColors(page, el("box plot viewer"), 2));
      await session.step(18, "When user sets \"Whisker Color\" property of box plot viewer to \"#1F77B4\"", () => setProperty(page, "Whisker Color", el("box plot viewer"), "#1F77B4"));
      await session.step(19, "Then \"Whisker Color\" property of box plot viewer should be \"#1F77B4\"", () => propertyShouldBe(page, "Whisker Color", el("box plot viewer"), "#1F77B4"));
      await session.step(20, "And box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(21, "When user sets \"Whisker Color\" property of box plot viewer to \"\"", () => setProperty(page, "Whisker Color", el("box plot viewer"), ""));
      await session.step(22, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
    });
    await run.scenario("The statistics strip and its ladder", async () => {
      await session.step(25, "Then \"Show Statistics\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Statistics", el("box plot viewer"), "true"));
      await session.step(26, "And box plot viewer should have a \"stats\" area", () => hasArea(page, el("box plot viewer"), "stats"));
      await session.step(27, "And close-stats icon in box plot viewer should be present", () => shouldBe(page, el("close-stats icon in box plot viewer"), "present"));
      await session.step(28, "And \"Show Avg\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Avg", el("box plot viewer"), "true"));
      await session.step(29, "And \"Show Med\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Med", el("box plot viewer"), "true"));
      await session.step(30, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]));
      await session.step(37, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(38, "And properties of box plot viewer should be:", () => propertiesShouldBe(page, el("box plot viewer"), [["Show Total Count","true"],["Show Inliers Count","true"],["Show Outliers Count","true"],["Show Stdev","true"],["Show Q1","true"],["Show Q3","true"]]));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
      await session.step(46, "When user sets \"Statistics Format\" property of box plot viewer to \"#,##0.00\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "#,##0.00"));
      await session.step(47, "Then \"Statistics Format\" property of box plot viewer should be \"#,##0.00\"", () => propertyShouldBe(page, "Statistics Format", el("box plot viewer"), "#,##0.00"));
      await session.step(48, "And \"Show P Value\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "true"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
      await session.step(50, "When user sets \"Statistics Format\" property of box plot viewer to \"auto\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "auto"));
      await session.step(51, "And user picks \"Show Total Count\" from the context menu of the \"stats\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Show Total Count", "stats", el("box plot viewer")));
      await session.step(52, "Then \"Show Total Count\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Total Count", el("box plot viewer"), "false"));
      await session.step(53, "When user picks \"Show Total Count\" from the context menu of the \"stats\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Show Total Count", "stats", el("box plot viewer")));
      await session.step(54, "Then \"Show Total Count\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Total Count", el("box plot viewer"), "true"));
      await session.step(55, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Total Count","false"],["Show Inliers Count","false"],["Show Outliers Count","false"],["Show Stdev","false"],["Show Q1","false"],["Show Q3","false"]]));
    });
    await run.scenario("The T key toggles the p-value", async () => {
      await session.step(64, "When user sets \"Show P Value\" property of box plot viewer to \"false\"", () => setProperty(page, "Show P Value", el("box plot viewer"), "false"));
      await session.step(65, "And user clicks on empty plot space of box plot viewer", () => clickEmptySpace(page, el("box plot viewer")));
      await session.step(66, "And user presses t", () => pressKey(page, "t"));
      await session.step(67, "Then \"Show P Value\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "true"));
      await session.step(68, "When user presses t", () => pressKey(page, "t"));
      await session.step(69, "Then \"Show P Value\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "false"));
      await session.step(70, "When user sets \"Show P Value\" property of box plot viewer to \"true\"", () => setProperty(page, "Show P Value", el("box plot viewer"), "true"));
    });
    await run.scenario("Three groups take the Alexander-Govern branch", async () => {
      await session.step(73, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(74, "And user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(75, "Then box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
      await session.step(76, "And box plot viewer should have a \"category Caucasian\" area", () => hasArea(page, el("box plot viewer"), "category Caucasian"));
      await session.step(77, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(78, "Then show-group-stats icon in box plot viewer should be visible", () => shouldBe(page, el("show-group-stats icon in box plot viewer"), "visible"));
      await session.step(79, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(80, "Then show-group-stats icon in box plot viewer should be hidden", () => shouldBe(page, el("show-group-stats icon in box plot viewer"), "hidden"));
      await session.step(81, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","SEX"],["Marker Color Column",""]]));
    });
    await run.scenario("The violin style", async () => {
      await session.step(86, "When user sets \"Plot Style\" property of box plot viewer to \"violin\"", () => setProperty(page, "Plot Style", el("box plot viewer"), "violin"));
      await session.step(87, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(88, "And the \"M values\" area of box plot viewer should be painted", () => areaPainted(page, "M values", el("box plot viewer")));
      await session.step(89, "And the \"F values\" area of box plot viewer should be painted", () => areaPainted(page, "F values", el("box plot viewer")));
      await session.step(90, "When user sets \"Bins\" property of box plot viewer to \"50\"", () => setProperty(page, "Bins", el("box plot viewer"), "50"));
      await session.step(91, "And user sets \"Bins\" property of box plot viewer to \"500\"", () => setProperty(page, "Bins", el("box plot viewer"), "500"));
      await session.step(92, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(93, "When user sets \"Interquartile Line Width\" property of box plot viewer to \"10\"", () => setProperty(page, "Interquartile Line Width", el("box plot viewer"), "10"));
      await session.step(94, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(95, "When user sets \"Violin Line Width\" property of box plot viewer to \"4\"", () => setProperty(page, "Violin Line Width", el("box plot viewer"), "4"));
      await session.step(96, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(97, "When user sets \"Violin Whisker Color\" property of box plot viewer to \"#00AA00\"", () => setProperty(page, "Violin Whisker Color", el("box plot viewer"), "#00AA00"));
      await session.step(98, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(99, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Plot Style","box"],["Bins","100"],["Interquartile Line Width","6"],["Violin Line Width","2"]]));
      await session.step(104, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
    });
    await run.scenario("Column color coding drives the marker colors", async () => {
      await session.step(107, "When user colors \"WEIGHT\" column linearly from \"#0000FF\" to \"#FF0000\"", () => colorLinear(page, "WEIGHT", "#0000FF", "#FF0000"));
      await session.step(108, "And user sets \"Marker Color Column\" property of box plot viewer to \"WEIGHT\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "WEIGHT"));
      await session.step(109, "Then \"Marker Color Column\" property of box plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "WEIGHT"));
      await session.step(110, "When user colors \"WEIGHT\" column linearly from \"#0000FF\" to \"#FF0000\" over 60 to 120", () => colorLinearOver(page, "WEIGHT", "#0000FF", "#FF0000", 60, 120));
      await session.step(111, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(112, "When user clears the row selection", () => clearSelection(page));
      await session.step(113, "And user drags a selection box over the \"M values\" area of box plot viewer", () => dragSelectionOverArea(page, "M values", el("box plot viewer")));
      await session.step(114, "Then some rows should be selected", () => someSelected(page));
      await session.step(115, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
      await session.step(116, "When user colors \"WEIGHT\" column conditionally:", () => colorConditional(page, "WEIGHT", [["50-90","#00FF00"],["90-150","#800080"]]));
      await session.step(119, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(120, "And box plot viewer should show a selection highlight", () => someHighlight(page, el("box plot viewer")));
      await session.step(121, "When user clears the row selection", () => clearSelection(page));
      await session.step(122, "And user removes the coloring of \"WEIGHT\" column", () => colorOff(page, "WEIGHT"));
      await session.step(123, "And user sets \"Marker Color Column\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "SEX"));
      await session.step(124, "And user colors \"SEX\" column categorically:", () => colorCategorical(page, "SEX", [["M","#E41A1C"]]));
      await session.step(126, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(127, "When user removes the coloring of \"SEX\" column", () => colorOff(page, "SEX"));
      await session.step(128, "And user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
    });
    await run.scenario("A datetime value", async () => {
      await session.step(131, "When user sets \"Value\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Value", el("box plot viewer"), "STARTED"));
      await session.step(132, "Then \"Value\" property of box plot viewer should be \"STARTED\"", () => propertyShouldBe(page, "Value", el("box plot viewer"), "STARTED"));
      await session.step(133, "And box plot viewer should be painted", () => painted(page, el("box plot viewer")));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
      await session.step(135, "When user sets \"Statistics Format\" property of box plot viewer to \"yyyy-MM-dd HH:mm\"", () => setProperty(page, "Statistics Format", el("box plot viewer"), "yyyy-MM-dd HH:mm"));
      await session.step(136, "Then \"Statistics Format\" property of box plot viewer should be \"yyyy-MM-dd HH:mm\"", () => propertyShouldBe(page, "Statistics Format", el("box plot viewer"), "yyyy-MM-dd HH:mm"));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
      await session.step(138, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","AGE"],["Statistics Format","auto"]]));
    });
    run.finish();
  });
});
