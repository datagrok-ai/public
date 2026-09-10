/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-density.feature
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaLessInk, areaMoreInk, areaPainted, areaRepainted, hasArea, hasNoArea, lessInk, moreInk, noErrors, propertyShouldBe, readingIs, readingReads, repainted, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot density overlay", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-density.feature", import.meta.url);
  test("PC plot density overlay", {tag: ["@journey", "@viewers", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"],["Show All Lines","false"]]));
    await session.step(22, "Then pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(23, "And the \"lines drawn\" reading of pc plot viewer should be 1", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1));
    await session.step(24, "And \"Show Density\" property of pc plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Density", el("pc plot viewer"), "false"));
    await session.step(25, "And the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
    await session.step(26, "And pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
    await run.scenario("The density style defaults to circles and Show Density draws them", async () => {
      await session.step(29, "Then \"Density Style\" property of pc plot viewer should be \"circles\"", () => propertyShouldBe(page, "Density Style", el("pc plot viewer"), "circles"));
      await session.step(30, "When user sets \"Show Density\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "true"));
      await session.step(31, "Then the \"density style\" reading of pc plot viewer should be \"circles\"", () => readingReads(page, "density style", el("pc plot viewer"), "circles"));
      await session.step(32, "And pc plot viewer should have a \"density \\\"AGE\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(33, "And pc plot viewer should have a \"density \\\"HEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"HEIGHT\""));
      await session.step(34, "And pc plot viewer should have a \"density \\\"WEIGHT\\\"\" area", () => hasArea(page, el("pc plot viewer"), "density \"WEIGHT\""));
      await session.step(35, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(36, "And pc plot viewer should have more ink than before", () => moreInk(page, el("pc plot viewer")));
      await session.step(37, "When user sets \"Show Density\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "false"));
      await session.step(38, "Then pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(39, "And pc plot viewer should have less ink than before", () => lessInk(page, el("pc plot viewer")));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Each density style draws a different shape in the same box", async () => {
      await session.step(43, "When user sets \"Show Density\" property of pc plot viewer to \"true\"", () => setProperty(page, "Show Density", el("pc plot viewer"), "true"));
      await session.step(44, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(45, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(46, "Then the \"density style\" reading of pc plot viewer should be \"box plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "box plot"));
      await session.step(47, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(48, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(49, "When user sets \"Density Style\" property of pc plot viewer to \"violin plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "violin plot"));
      await session.step(50, "Then the \"density style\" reading of pc plot viewer should be \"violin plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "violin plot"));
      await session.step(51, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(52, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(53, "When user sets \"Density Style\" property of pc plot viewer to \"circles\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "circles"));
      await session.step(54, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The box plot's parts each add their own ink", async () => {
      await session.step(58, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(59, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(60, "When user sets \"Show Circles\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Circles", el("pc plot viewer"), "false"));
      await session.step(61, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(62, "When user sets \"Show Interquartile Range\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Interquartile Range", el("pc plot viewer"), "false"));
      await session.step(63, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(64, "When user sets \"Show Median\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Median", el("pc plot viewer"), "false"));
      await session.step(65, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(66, "When user sets \"Show Mean Cross\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Mean Cross", el("pc plot viewer"), "false"));
      await session.step(67, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(68, "When user sets \"Show Upper Dash\" property of pc plot viewer to \"false\"", () => setProperty(page, "Show Upper Dash", el("pc plot viewer"), "false"));
      await session.step(69, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have less ink than before", () => areaLessInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(70, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Upper Dash","true"],["Show Mean Cross","true"],["Show Median","true"],["Show Interquartile Range","true"],["Show Circles","true"]]));
      await session.step(76, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have more ink than before", () => areaMoreInk(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bin count changes the violin's outline", async () => {
      await session.step(80, "When user sets \"Density Style\" property of pc plot viewer to \"violin plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "violin plot"));
      await session.step(81, "Then \"Bins\" property of pc plot viewer should be \"100\"", () => propertyShouldBe(page, "Bins", el("pc plot viewer"), "100"));
      await session.step(82, "When user sets \"Bins\" property of pc plot viewer to \"200\"", () => setProperty(page, "Bins", el("pc plot viewer"), "200"));
      await session.step(83, "Then \"Bins\" property of pc plot viewer should be \"200\"", () => propertyShouldBe(page, "Bins", el("pc plot viewer"), "200"));
      await session.step(84, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(85, "When user sets \"Bins\" property of pc plot viewer to \"20\"", () => setProperty(page, "Bins", el("pc plot viewer"), "20"));
      await session.step(86, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(87, "When user sets \"Bins\" property of pc plot viewer to \"100\"", () => setProperty(page, "Bins", el("pc plot viewer"), "100"));
      await session.step(88, "Then the \"density \\\"AGE\\\"\" area of pc plot viewer should have repainted", () => areaRepainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The density survives a normalization double-toggle (github-1546)", async () => {
      await session.step(92, "When user sets \"Density Style\" property of pc plot viewer to \"box plot\"", () => setProperty(page, "Density Style", el("pc plot viewer"), "box plot"));
      await session.step(93, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"false\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "false"));
      await session.step(94, "Then the \"normalization\" reading of pc plot viewer should be \"global\"", () => readingReads(page, "normalization", el("pc plot viewer"), "global"));
      await session.step(95, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(96, "When user sets \"Normalize Each Column\" property of pc plot viewer to \"true\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "true"));
      await session.step(97, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"false\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "false"));
      await session.step(98, "And user sets \"Normalize Each Column\" property of pc plot viewer to \"true\"", () => setProperty(page, "Normalize Each Column", el("pc plot viewer"), "true"));
      await session.step(99, "Then the \"normalization\" reading of pc plot viewer should be \"per column\"", () => readingReads(page, "normalization", el("pc plot viewer"), "per column"));
      await session.step(100, "And the \"density style\" reading of pc plot viewer should be \"box plot\"", () => readingReads(page, "density style", el("pc plot viewer"), "box plot"));
      await session.step(101, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(102, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic AGE axis keeps the density painted", async () => {
      await session.step(106, "When user sets \"Log Columns\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Log Columns", el("pc plot viewer"), "AGE"));
      await session.step(107, "Then \"Log Columns\" property of pc plot viewer should be \"AGE\"", () => propertyShouldBe(page, "Log Columns", el("pc plot viewer"), "AGE"));
      await session.step(108, "And pc plot viewer should have repainted", () => repainted(page, el("pc plot viewer")));
      await session.step(109, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(110, "And the \"density \\\"WEIGHT\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"WEIGHT\"", el("pc plot viewer")));
      await session.step(111, "When user sets \"Log Columns\" property of pc plot viewer to \"\"", () => setProperty(page, "Log Columns", el("pc plot viewer"), ""));
      await session.step(112, "Then \"Log Columns\" property of pc plot viewer should be \"\"", () => propertyShouldBe(page, "Log Columns", el("pc plot viewer"), ""));
      await session.step(113, "And the \"density \\\"AGE\\\"\" area of pc plot viewer should be painted", () => areaPainted(page, "density \"AGE\"", el("pc plot viewer")));
      await session.step(114, "And pc plot viewer should have repainted", () => repainted(page, el("pc plot viewer")));
      await session.step(115, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Density","false"],["Density Style","circles"],["Show All Lines","true"]]));
      await session.step(119, "Then the \"density style\" reading of pc plot viewer should be \"\"", () => readingReads(page, "density style", el("pc plot viewer"), ""));
      await session.step(120, "And pc plot viewer should not have a \"density \\\"AGE\\\"\" area", () => hasNoArea(page, el("pc plot viewer"), "density \"AGE\""));
      await session.step(121, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
