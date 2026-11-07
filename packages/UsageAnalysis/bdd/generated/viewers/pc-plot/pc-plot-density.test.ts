/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-density.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot, core#967]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, applyNamedLayout, areaLessInk, areaMoreInk, areaPainted, areaRepainted, hasArea, hasNoArea, hoverArea, lessInk, moreInk, noErrors, painted, pointerAway, propertyShouldBe, readingIs, readingReads, repainted, saveNamedLayout, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot density overlay", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-density.feature", import.meta.url);
  test("PC plot density overlay", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot", "@edge-cases", "@realizes:core#967"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(22, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"],["Show All Lines","false"]]), [["Column Names","AGE, HEIGHT, WEIGHT"],["Show All Lines","false"]]);
    await session.step(25, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(26, "And the \"lines drawn\" reading of pc plot viewer should be 1", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1));
    await session.step(27, "And \"Show Density\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Density", el("pc plot viewer"), "false"));
    await session.step(28, "And the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
    await session.step(29, "And pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
    await run.scenario("The density style defaults to circles and Show Density draws them", async () => {
      await session.step(32, "Then \"Density Style\" property of pc plot viewer should be \"circles\"", () => propertyShouldBe(page, "Density Style", el("pc plot viewer"), "circles"));
      await session.step(33, "When user sets \"Show Density\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "true"));
      await session.step(34, "Then the \"density style\" reading of pc plot viewer should be \"circles\"", () => readingReads(page, "density style", el("pc plot viewer"), "circles"));
      await session.step(35, "And pc plot viewer should have a \"density \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(36, "And pc plot viewer should have a \"density \\\"HEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"HEIGHT\""));
      await session.step(37, "And pc plot viewer should have a \"density \\\"WEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"WEIGHT\""));
      await session.step(38, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(39, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(40, "When user sets \"Show Density\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "false"));
      await session.step(41, "Then pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(42, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Each density style draws a different shape in the same box", async () => {
      await session.step(46, "When user sets \"Show Density\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "true"));
      await session.step(47, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(48, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(49, "Then the \"density style\" reading of pc plot viewer should be \"box plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "box plot"));
      await session.step(50, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(51, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(52, "When user sets \"Density Style\" property of pc plot viewer to \"violin plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "violin plot"));
      await session.step(53, "Then the \"density style\" reading of pc plot viewer should be \"violin plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "violin plot"));
      await session.step(54, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(55, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(56, "When user sets \"Density Style\" property of pc plot viewer to \"circles\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "circles"));
      await session.step(57, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The box plot's parts each add their own ink", async () => {
      await session.step(61, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(62, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(63, "When user sets \"Show Circles\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Circles", el("pc plot viewer"), "false"));
      await session.step(64, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(65, "When user sets \"Show Interquartile Range\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Interquartile Range", el("pc plot viewer"), "false"));
      await session.step(66, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(67, "When user sets \"Show Median\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Median", el("pc plot viewer"), "false"));
      await session.step(68, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(69, "When user sets \"Show Mean Cross\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Mean Cross", el("pc plot viewer"), "false"));
      await session.step(70, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(71, "When user sets \"Show Upper Dash\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Upper Dash", el("pc plot viewer"), "false"));
      await session.step(72, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(73, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Upper Dash","true"],["Show Mean Cross","true"],["Show Median","true"],["Show Interquartile Range","true"],["Show Circles","true"]]), [["Show Upper Dash","true"],["Show Mean Cross","true"],["Show Median","true"],["Show Interquartile Range","true"],["Show Circles","true"]]);
      await session.step(79, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have more ink than before", () => areaMoreInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bin count changes the violin's outline", async () => {
      await session.step(83, "When user sets \"Density Style\" property of pc plot viewer to \"violin plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "violin plot"));
      await session.step(84, "Then \"Bins\" property of pc plot viewer should be \"100\"", () => propertyShouldBe(page, "Bins", el("pc plot viewer"), "100"));
      await session.step(85, "When user sets \"Bins\" property of pc plot viewer to \"200\"", () => setProperty(page, "Bins", el("pc plot viewer"), "200"));
      await session.step(86, "Then \"Bins\" property of pc plot viewer should be \"200\"", () => propertyShouldBe(page, "Bins", el("pc plot viewer"), "200"));
      await session.step(87, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(88, "When user sets \"Bins\" property of pc plot viewer to \"20\"", () => setProperty(page, "Bins", el("pc plot viewer"), "20"));
      await session.step(89, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(90, "When user sets \"Bins\" property of pc plot viewer to \"100\"", () => setProperty(page, "Bins", el("pc plot viewer"), "100"));
      await session.step(91, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The density survives a normalization double-toggle (github-1546)", async () => {
      await session.step(95, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(96, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"false\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "false"));
      await session.step(97, "Then the \"normalization\" reading of pc plot viewer should be \"global\"", () => readingReads(page, "normalization", el("pc plot viewer"), "global"));
      await session.step(98, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(99, "When user sets \"Normalize Each Column\" property of pc plot viewer to \"true\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "true"));
      await session.step(100, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"false\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "false"));
      await session.step(101, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"true\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "true"));
      await session.step(102, "Then the \"normalization\" reading of pc plot viewer should be \"per column\"", () => readingReads(page, "normalization", el("pc plot viewer"), "per column"));
      await session.step(103, "And the \"density style\" reading of pc plot viewer should be \"box plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "box plot"));
      await session.step(104, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(105, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
      await session.step(106, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic AGE axis keeps the density painted", async () => {
      await session.step(109, "When user sets \"Log Columns\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Log Columns", el("pc plot viewer"), "AGE"));
      await session.step(110, "Then \"Log Columns\" property of pc plot viewer should be \"AGE\"", () => propertyShouldBe(page, "Log Columns", el("pc plot viewer"), "AGE"));
      await session.step(111, "And pc plot viewer should have repainted", () => repainted(page, el("pc plot viewer")));
      await session.step(112, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(113, "And the \"density \\\"WEIGHT\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"WEIGHT\"", el("pc plot viewer")));
      await session.step(114, "When user sets \"Log Columns\" property of pc plot viewer to \"\"", () => setProperty(page, "Log Columns", el("pc plot viewer"), ""));
      await session.step(115, "Then \"Log Columns\" property of pc plot viewer should be \"\"", () => propertyShouldBe(page, "Log Columns", el("pc plot viewer"), ""));
      await session.step(116, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(117, "And pc plot viewer should have repainted", () => repainted(page, el("pc plot viewer")));
      await session.step(118, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Density","false"],["Density Style","circles"],["Show All Lines","true"]]), [["Show Density","false"],["Density Style","circles"],["Show All Lines","true"]]);
      await session.step(122, "Then the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
      await session.step(123, "And pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(124, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved with Line Width 0 and density on is applied and hovered without an error (core#967)", async () => {
      await session.step(129, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Density","true"],["Line Width","0"]]), [["Show Density","true"],["Line Width","0"]]);
      await session.step(132, "And user saves the layout of the current table view as \"pc-width-0\"", () => saveNamedLayout(page, "pc-width-0"));
      await session.step(133, "And user applies the layout \"pc-width-0\" to the current table view", () => applyNamedLayout(page, "pc-width-0"));
      await session.step(134, "Then \"Line Width\" property of pc plot viewer should be \"0\"", () => propertyShouldBe(page, "Line Width", el("pc plot viewer"), "0"));
      await session.step(135, "And the \"density style\" reading of pc plot viewer should be \"circles\"", () => readingReads(page, "density style", el("pc plot viewer"), "circles"));
      await session.step(136, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
      await session.step(137, "When user hovers over the \"view\" area of pc plot viewer", () => hoverArea(page, "view", el("pc plot viewer")));
      await session.step(138, "And user moves the pointer away from pc plot viewer", () => pointerAway(page, el("pc plot viewer")));
      await session.step(139, "Then the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
      await session.step(141, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Density","false"],["Line Width","1"]]), [["Show Density","false"],["Line Width","1"]]);
      await session.step(144, "Then the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
      await session.step(145, "And pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
