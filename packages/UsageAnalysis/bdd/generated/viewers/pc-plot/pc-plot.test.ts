/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {axesShouldBe, dragAxisLabel} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hoverOver, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaNarrower, areaShorter, areaTaller, areaWider, hasArea, hasNoArea, lessInk, moreInk, noErrors, painted, pickFromContextMenu, propertyShouldBe, readingIs, readingReads, repaintedBy, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot axes, vertical scale and chrome", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot.feature", import.meta.url);
  test("PC plot axes, vertical scale and chrome", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(25, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
    await session.step(27, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(28, "And the \"axes\" reading of pc plot viewer should be 3", () => readingIs(page, "axes", el("pc plot viewer"), 3));
    await session.step(29, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(30, "And the \"normalization\" reading of pc plot viewer should be \"per column\"", () => readingReads(page, "normalization", el("pc plot viewer"), "per column"));
    await session.step(31, "And the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
    await session.step(32, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(33, "And pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
    await run.scenario("One axis per column, in the order the plot draws them", async () => {
      await session.step(36, "Then pc plot viewer should have an \"axis \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis \"AGE\""));
      await session.step(37, "And pc plot viewer should have an \"axis \\\"HEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis \"HEIGHT\""));
      await session.step(38, "And pc plot viewer should have an \"axis \\\"WEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis \"WEIGHT\""));
      await session.step(39, "And pc plot viewer should have a \"band \\\"AGE\\\" - \\\"HEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"AGE\" - \"HEIGHT\""));
      await session.step(40, "And pc plot viewer should have a \"band \\\"HEIGHT\\\" - \\\"WEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"HEIGHT\" - \"WEIGHT\""));
      await session.step(41, "And pc plot viewer should not have a \"band \\\"AGE\\\" - \\\"WEIGHT\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "band \"AGE\" - \"WEIGHT\""));
      await session.step(42, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(43, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Adding a column adds its axis and removing it takes it back (GROK-18000)", async () => {
      await session.step(47, "When user sets \"Column Names\" property of pc plot viewer to \"AGE, HEIGHT, WEIGHT, STARTED\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(48, "Then the \"axes\" reading of pc plot viewer should be 4", () => readingIs(page, "axes", el("pc plot viewer"), 4));
      await session.step(49, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(50, "And pc plot viewer should have an \"axis \\\"STARTED\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis \"STARTED\""));
      await session.step(51, "And pc plot viewer should have a \"band \\\"WEIGHT\\\" - \\\"STARTED\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"WEIGHT\" - \"STARTED\""));
      await session.step(52, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(53, "When user sets \"Column Names\" property of pc plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(54, "Then the \"axes\" reading of pc plot viewer should be 3", () => readingIs(page, "axes", el("pc plot viewer"), 3));
      await session.step(55, "And the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(56, "And pc plot viewer should not have an \"axis \\\"STARTED\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis \"STARTED\""));
      await session.step(57, "And pc plot viewer should not have a \"band \\\"WEIGHT\\\" - \\\"STARTED\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "band \"WEIGHT\" - \"STARTED\""));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Y Axis > Global draws every axis to one range, Normalized to its column's own", async () => {
      await session.step(61, "Then the \"axis min of \\\"AGE\\\"\" reading of pc plot viewer should be 18", () => readingIs(page, "axis min of \"AGE\"", el("pc plot viewer"), 18));
      await session.step(62, "And the \"axis max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "axis max of \"AGE\"", el("pc plot viewer"), 89));
      await session.step(63, "And the \"axis min of \\\"HEIGHT\\\"\" reading of pc plot viewer should be 137.32200622558594", () => readingIs(page, "axis min of \"HEIGHT\"", el("pc plot viewer"), 137.32200622558594));
      await session.step(64, "And the \"global min\" reading of pc plot viewer should be 18", () => readingIs(page, "global min", el("pc plot viewer"), 18));
      await session.step(65, "And the \"global max\" reading of pc plot viewer should be 198.86199951171875", () => readingIs(page, "global max", el("pc plot viewer"), 198.86199951171875));
      await session.step(66, "And pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(67, "And pc plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("pc plot viewer"), "y axis"));
      await session.step(68, "When user picks \"Y Axis > Global\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Y Axis > Global", el("pc plot viewer")));
      await session.step(69, "Then \"Normalize Each Column\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Normalize Each Column", el("pc plot viewer"), "false"));
      await session.step(70, "And the \"normalization\" reading of pc plot viewer should be \"global\"", () => readingReads(page, "normalization", el("pc plot viewer"), "global"));
      await session.step(71, "And the \"axis max of \\\"AGE\\\"\" reading of pc plot viewer should be 198.86199951171875", () => readingIs(page, "axis max of \"AGE\"", el("pc plot viewer"), 198.86199951171875));
      await session.step(72, "And the \"axis min of \\\"HEIGHT\\\"\" reading of pc plot viewer should be 18", () => readingIs(page, "axis min of \"HEIGHT\"", el("pc plot viewer"), 18));
      await session.step(73, "And pc plot viewer should have a \"y axis\" area", () => hasArea(page, el("pc plot viewer"), "y axis"));
      await session.step(74, "And pc plot viewer should not have an \"axis min \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(75, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(76, "When user picks \"Y Axis > Normalized\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Y Axis > Normalized", el("pc plot viewer")));
      await session.step(77, "Then \"Normalize Each Column\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Normalize Each Column", el("pc plot viewer"), "true"));
      await session.step(78, "And the \"normalization\" reading of pc plot viewer should be \"per column\"", () => readingReads(page, "normalization", el("pc plot viewer"), "per column"));
      await session.step(79, "And the \"axis max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "axis max of \"AGE\"", el("pc plot viewer"), 89));
      await session.step(80, "And the \"axis min of \\\"HEIGHT\\\"\" reading of pc plot viewer should be 137.32200622558594", () => readingIs(page, "axis min of \"HEIGHT\"", el("pc plot viewer"), 137.32200622558594));
      await session.step(81, "And pc plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("pc plot viewer"), "y axis"));
      await session.step(82, "And pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(83, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging a column label reorders the axes", async () => {
      await session.step(87, "When user drags the \"WEIGHT\" axis label of pc plot viewer onto the \"AGE\" axis label", () => dragAxisLabel(page, "WEIGHT", "AGE"));
      await session.step(88, "Then the axes of pc plot viewer should be \"WEIGHT, AGE, HEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "WEIGHT, AGE, HEIGHT"));
      await session.step(89, "And \"Column Names\" property of pc plot viewer should be \"WEIGHT, AGE, HEIGHT\"", () => propertyShouldBe(page, "Column Names", el("pc plot viewer"), "WEIGHT, AGE, HEIGHT"));
      await session.step(90, "And pc plot viewer should have a \"band \\\"WEIGHT\\\" - \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "band \"WEIGHT\" - \"AGE\""));
      await session.step(91, "And pc plot viewer should not have a \"band \\\"HEIGHT\\\" - \\\"WEIGHT\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "band \"HEIGHT\" - \"WEIGHT\""));
      await session.step(92, "When user sets \"Column Names\" property of pc plot viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(93, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(94, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Title and description show and clear", async () => {
      await session.step(97, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Title","true"],["Title","Demographics"]]));
      await session.step(100, "Then title of pc plot viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of pc plot viewer"), "Demographics"));
      await session.step(101, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Description","Three axes"],["Description Visibility Mode","Always"]]));
      await session.step(104, "Then description of pc plot viewer should have text \"Three axes\"", () => shouldHaveText(page, el("description of pc plot viewer"), "Three axes"));
      await session.step(105, "When user sets \"Description Position\" property of pc plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("pc plot viewer"), "Bottom"));
      await session.step(106, "Then description of pc plot viewer should be visible", () => shouldBe(page, el("description of pc plot viewer"), "visible"));
      await session.step(107, "When user sets \"Description Visibility Mode\" property of pc plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("pc plot viewer"), "Never"));
      await session.step(108, "Then description of pc plot viewer should be absent", () => shouldBe(page, el("description of pc plot viewer"), "absent"));
      await session.step(109, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
      await session.step(115, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Line width, label orientation and the horizontal margin change what is drawn", async () => {
      await session.step(118, "When user sets \"Line Width\" property of pc plot viewer to \"3\"", () => setProperty(page, "Line Width", el("pc plot viewer"), "3"));
      await session.step(119, "Then pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(120, "When user sets \"Line Width\" property of pc plot viewer to \"1\"", () => setProperty(page, "Line Width", el("pc plot viewer"), "1"));
      await session.step(121, "Then pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(122, "When user sets \"Current Line Width\" property of pc plot viewer to \"8\"", () => setProperty(page, "Current Line Width", el("pc plot viewer"), "8"));
      await session.step(123, "Then the \"current row\" reading of pc plot viewer should be 1", () => readingIs(page, "current row", el("pc plot viewer"), 1));
      await session.step(124, "And pc plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pc plot viewer"), 500));
      await session.step(125, "When user sets \"Current Line Width\" property of pc plot viewer to \"2\"", () => setProperty(page, "Current Line Width", el("pc plot viewer"), "2"));
      await session.step(126, "Then pc plot viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pc plot viewer"), 500));
      await session.step(127, "When user sets \"Horz Margin\" property of pc plot viewer to \"90\"", () => setProperty(page, "Horz Margin", el("pc plot viewer"), "90"));
      await session.step(128, "Then the \"view\" area of pc plot viewer should be narrower than before", () => areaNarrower(page, "view", el("pc plot viewer")));
      await session.step(129, "And pc plot viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pc plot viewer"), 1000));
      await session.step(130, "When user sets \"Horz Margin\" property of pc plot viewer to \"40\"", () => setProperty(page, "Horz Margin", el("pc plot viewer"), "40"));
      await session.step(131, "Then the \"view\" area of pc plot viewer should be wider than before", () => areaWider(page, "view", el("pc plot viewer")));
      await session.step(132, "When user sets \"Labels Orientation\" property of pc plot viewer to \"Vert\"", () => setProperty(page, "Labels Orientation", el("pc plot viewer"), "Vert"));
      await session.step(133, "Then the \"view\" area of pc plot viewer should be shorter than before", () => areaShorter(page, "view", el("pc plot viewer")));
      await session.step(134, "When user sets \"Labels Orientation\" property of pc plot viewer to \"Auto\"", () => setProperty(page, "Labels Orientation", el("pc plot viewer"), "Auto"));
      await session.step(135, "Then the \"view\" area of pc plot viewer should be taller than before", () => areaTaller(page, "view", el("pc plot viewer")));
      await session.step(136, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Min Max and Show Labels drop the axis chrome", async () => {
      await session.step(139, "When user sets \"Show Min Max\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Min Max", el("pc plot viewer"), "false"));
      await session.step(140, "Then pc plot viewer should not have an \"axis min \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(141, "And pc plot viewer should not have an \"axis max \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis max \"AGE\""));
      await session.step(142, "And pc plot viewer should have repainted by at least 200 pixels", () => repaintedBy(page, el("pc plot viewer"), 200));
      await session.step(143, "When user sets \"Show Labels\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Labels", el("pc plot viewer"), "false"));
      await session.step(144, "Then pc plot viewer should not have an \"axis label \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(145, "And pc plot viewer should have repainted by at least 200 pixels", () => repaintedBy(page, el("pc plot viewer"), 200));
      await session.step(146, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Min Max","true"],["Show Labels","true"]]));
      await session.step(149, "Then pc plot viewer should have an \"axis min \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis min \"AGE\""));
      await session.step(150, "And pc plot viewer should have an \"axis label \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "axis label \"AGE\""));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Filter > Show Filters hides the in-chart range sliders", async () => {
      await session.step(154, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(155, "Then pc plot viewer should have a \"range slider \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range slider \"AGE\""));
      await session.step(156, "And pc plot viewer should have a \"range max handle \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range max handle \"AGE\""));
      await session.step(157, "And pc plot viewer should have a \"range slider \\\"WEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range slider \"WEIGHT\""));
      await session.step(158, "When user picks \"Filter > Show Filters\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filters", el("pc plot viewer")));
      await session.step(159, "Then \"Show Filters\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Filters", el("pc plot viewer"), "false"));
      await session.step(160, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(161, "Then pc plot viewer should not have a \"range slider \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "range slider \"AGE\""));
      await session.step(162, "And pc plot viewer should not have a \"range max handle \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "range max handle \"AGE\""));
      await session.step(163, "And the \"range max of \\\"AGE\\\"\" reading of pc plot viewer should be 89", () => readingIs(page, "range max of \"AGE\"", el("pc plot viewer"), 89));
      await session.step(164, "When user picks \"Filter > Show Filters\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "Filter > Show Filters", el("pc plot viewer")));
      await session.step(165, "Then \"Show Filters\" property of pc plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Filters", el("pc plot viewer"), "true"));
      await session.step(166, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
      await session.step(167, "Then pc plot viewer should have a \"range slider \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "range slider \"AGE\""));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("To Script > To JavaScript prints the call that rebuilds the plot", async () => {
      await session.step(171, "Given the package autostarts have completed", () => autostartsCompleted(page));
      await session.step(172, "When user picks \"To Script > To JavaScript\" from the context menu of pc plot viewer", () => pickFromContextMenu(page, "To Script > To JavaScript", el("pc plot viewer")));
      await session.step(173, "Then balloon should contain text \"addViewer\"", () => shouldContainText(page, el("balloon"), "addViewer"));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
