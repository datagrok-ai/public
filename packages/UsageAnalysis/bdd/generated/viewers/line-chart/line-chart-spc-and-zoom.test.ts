/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-spc-and-zoom.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaLessInk, areaPainted, areaShorter, eventFired, hasArea, hasNoArea, listenFor, noErrors, pickFromContextMenu, propertyShouldBe, readingBetween, readingHigher, readingIs, readingLower, repainted, reportsNoError, setProperties, setProperty, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart statistical process control, zoom and Reset View", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-spc-and-zoom.feature", import.meta.url);
  test("Line chart statistical process control, zoom and Reset View", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(22, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]));
    await session.step(25, "Then 100 rows should pass the filter", () => filterPasses(page, 100));
    await session.step(26, "And \"showStatisticalProcessControl\" property of line chart viewer should be \"false\"", () => propertyShouldBe(page, "showStatisticalProcessControl", el("line chart viewer"), "false"));
    await session.step(27, "And line chart viewer should not have a \"control limits\" area", () => hasNoArea(page, el("line chart viewer"), "control limits"));
    await session.step(28, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
    await session.step(29, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Switching SPC on draws the limits, the sigma bands and the average", async () => {
      await session.step(32, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"true\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "true"));
      await session.step(33, "Then line chart viewer should have a \"control limits\" area", () => hasArea(page, el("line chart viewer"), "control limits"));
      await session.step(34, "And line chart viewer should have a \"sigma 1\" area", () => hasArea(page, el("line chart viewer"), "sigma 1"));
      await session.step(35, "And line chart viewer should have a \"sigma 2\" area", () => hasArea(page, el("line chart viewer"), "sigma 2"));
      await session.step(36, "And line chart viewer should have an \"average\" area", () => hasArea(page, el("line chart viewer"), "average"));
      await session.step(37, "And the \"spc average\" reading of line chart viewer should be between 3.4 and 3.6", () => readingBetween(page, "spc average", el("line chart viewer"), 3.4, 3.6));
      await session.step(38, "And the \"upper control limit\" reading of line chart viewer should be between 21 and 22", () => readingBetween(page, "upper control limit", el("line chart viewer"), 21, 22));
      await session.step(39, "And the \"lower control limit\" reading of line chart viewer should be between -15 and -14", () => readingBetween(page, "lower control limit", el("line chart viewer"), -15, -14));
      await session.step(40, "And the \"violations\" reading of line chart viewer should be 0", () => readingIs(page, "violations", el("line chart viewer"), 0));
      await session.step(41, "And the \"control limits\" area of line chart viewer should be painted", () => areaPainted(page, "control limits", el("line chart viewer")));
      await session.step(42, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(43, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(44, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"false\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "false"));
      await session.step(45, "Then line chart viewer should not have a \"control limits\" area", () => hasNoArea(page, el("line chart viewer"), "control limits"));
      await session.step(46, "And line chart viewer should not have an \"average\" area", () => hasNoArea(page, el("line chart viewer"), "average"));
      await session.step(47, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be lower than before", () => readingLower(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(48, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("SPC is gated on a single un-split series", async () => {
      await session.step(52, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"true\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "true"));
      await session.step(53, "Then line chart viewer should have a \"control limits\" area", () => hasArea(page, el("line chart viewer"), "control limits"));
      await session.step(54, "When user sets \"splitColumnNames\" property of line chart viewer to \"Stereo Category\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), "Stereo Category"));
      await session.step(55, "Then line chart viewer should not have a \"control limits\" area", () => hasNoArea(page, el("line chart viewer"), "control limits"));
      await session.step(56, "And line chart viewer should not have a \"sigma 1\" area", () => hasNoArea(page, el("line chart viewer"), "sigma 1"));
      await session.step(57, "And line chart viewer should not have an \"average\" area", () => hasNoArea(page, el("line chart viewer"), "average"));
      await session.step(58, "And the \"lines\" reading of line chart viewer should be 5", () => readingIs(page, "lines", el("line chart viewer"), 5));
      await session.step(59, "When user sets \"splitColumnNames\" property of line chart viewer to \"\"", () => setProperty(page, "splitColumnNames", el("line chart viewer"), ""));
      await session.step(60, "Then line chart viewer should have a \"control limits\" area", () => hasArea(page, el("line chart viewer"), "control limits"));
      await session.step(61, "When user sets \"multiAxis\" property of line chart viewer to \"true\"", () => setProperty(page, "multiAxis", el("line chart viewer"), "true"));
      await session.step(62, "Then line chart viewer should not have a \"control limits\" area", () => hasNoArea(page, el("line chart viewer"), "control limits"));
      await session.step(63, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["multiAxis","false"],["showStatisticalProcessControl","false"]]));
      await session.step(66, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bands can be hidden one at a time and the limits stay reported", async () => {
      await session.step(69, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"true\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "true"));
      await session.step(70, "Then line chart viewer should have a \"sigma 1\" area", () => hasArea(page, el("line chart viewer"), "sigma 1"));
      await session.step(71, "When user sets \"showSigma1\" property of line chart viewer to \"false\"", () => setProperty(page, "showSigma1", el("line chart viewer"), "false"));
      await session.step(72, "Then line chart viewer should not have a \"sigma 1\" area", () => hasNoArea(page, el("line chart viewer"), "sigma 1"));
      await session.step(73, "And line chart viewer should have a \"sigma 2\" area", () => hasArea(page, el("line chart viewer"), "sigma 2"));
      await session.step(74, "And the \"plot\" area of line chart viewer should have less ink than before", () => areaLessInk(page, "plot", el("line chart viewer")));
      await session.step(75, "When user sets \"showSigma2\" property of line chart viewer to \"false\"", () => setProperty(page, "showSigma2", el("line chart viewer"), "false"));
      await session.step(76, "Then line chart viewer should not have a \"sigma 2\" area", () => hasNoArea(page, el("line chart viewer"), "sigma 2"));
      await session.step(77, "When user sets \"showAverage\" property of line chart viewer to \"false\"", () => setProperty(page, "showAverage", el("line chart viewer"), "false"));
      await session.step(78, "Then line chart viewer should not have an \"average\" area", () => hasNoArea(page, el("line chart viewer"), "average"));
      await session.step(79, "And the \"spc average\" reading of line chart viewer should be between 3.4 and 3.6", () => readingBetween(page, "spc average", el("line chart viewer"), 3.4, 3.6));
      await session.step(80, "When user sets \"showControlLimits\" property of line chart viewer to \"false\"", () => setProperty(page, "showControlLimits", el("line chart viewer"), "false"));
      await session.step(81, "Then line chart viewer should not have a \"control limits\" area", () => hasNoArea(page, el("line chart viewer"), "control limits"));
      await session.step(82, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be lower than before", () => readingLower(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(83, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showSigma1","true"],["showSigma2","true"],["showAverage","true"],["showControlLimits","true"]]));
      await session.step(88, "Then line chart viewer should have a \"control limits\" area", () => hasArea(page, el("line chart viewer"), "control limits"));
      await session.step(89, "And the 'y axis max of \"Chemical Space X\"' reading of line chart viewer should be higher than before", () => readingHigher(page, "y axis max of \"Chemical Space X\"", el("line chart viewer")));
      await session.step(90, "And the \"control limits\" area of line chart viewer should be painted", () => areaPainted(page, "control limits", el("line chart viewer")));
      await session.step(91, "When user sets \"showStatisticalProcessControl\" property of line chart viewer to \"false\"", () => setProperty(page, "showStatisticalProcessControl", el("line chart viewer"), "false"));
      await session.step(92, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Western Electric rules flag points, and hand-set limits flag far more of them", async () => {
      await session.step(95, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showStatisticalProcessControl","true"],["showBias","true"],["showConsistentTrend","true"],["showOscillation","true"],["showMediumShift","true"],["showSustainedShift","true"],["showSuppressedVariation","true"]]));
      await session.step(103, "Then the \"violations\" reading of line chart viewer should be 33", () => readingIs(page, "violations", el("line chart viewer"), 33));
      await session.step(104, "And line chart viewer should have a \"violation 2\" area", () => hasArea(page, el("line chart viewer"), "violation 2"));
      await session.step(105, "And line chart viewer should have a \"violation 6\" area", () => hasArea(page, el("line chart viewer"), "violation 6"));
      await session.step(106, "And line chart viewer should not have a \"violation 1\" area", () => hasNoArea(page, el("line chart viewer"), "violation 1"));
      await session.step(107, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["lowerControlLimit","0"],["upperControlLimit","5"]]));
      await session.step(110, "Then the \"upper control limit\" reading of line chart viewer should be 5", () => readingIs(page, "upper control limit", el("line chart viewer"), 5));
      await session.step(111, "And the \"lower control limit\" reading of line chart viewer should be 0", () => readingIs(page, "lower control limit", el("line chart viewer"), 0));
      await session.step(112, "And the \"violations\" reading of line chart viewer should be 89", () => readingIs(page, "violations", el("line chart viewer"), 89));
      await session.step(113, "And line chart viewer should have a \"violation 1\" area", () => hasArea(page, el("line chart viewer"), "violation 1"));
      await session.step(114, "And the \"control limits\" area of line chart viewer should be shorter than before", () => areaShorter(page, "control limits", el("line chart viewer")));
      await session.step(115, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["lowerControlLimit",""],["upperControlLimit",""]]));
      await session.step(118, "Then the \"upper control limit\" reading of line chart viewer should be between 21 and 22", () => readingBetween(page, "upper control limit", el("line chart viewer"), 21, 22));
      await session.step(119, "And the \"violations\" reading of line chart viewer should be 33", () => readingIs(page, "violations", el("line chart viewer"), 33));
      await session.step(120, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["showBias","false"],["showConsistentTrend","false"],["showOscillation","false"],["showMediumShift","false"],["showSustainedShift","false"],["showSuppressedVariation","false"],["showStatisticalProcessControl","false"]]));
      await session.step(128, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The wheel zooms the X axis and Reset View puts it back", async () => {
      await session.step(131, "Given user listens for \"d4-linechart-zoomed\" event on line chart viewer", () => listenFor(page, "d4-linechart-zoomed", el("line chart viewer")));
      await session.step(132, "And user listens for \"d4-linechart-reset-view\" event on line chart viewer", () => listenFor(page, "d4-linechart-reset-view", el("line chart viewer")));
      await session.step(133, "When user scrolls the mouse wheel up over the \"plot\" area of line chart viewer", () => wheelOverArea(page, "up", "plot", el("line chart viewer")));
      await session.step(134, "Then \"d4-linechart-zoomed\" event should have fired on line chart viewer", () => eventFired(page, "d4-linechart-zoomed", el("line chart viewer")));
      await session.step(135, "And the \"x axis span\" reading of line chart viewer should be lower than before", () => readingLower(page, "x axis span", el("line chart viewer")));
      await session.step(136, "And line chart viewer should have repainted", () => repainted(page, el("line chart viewer")));
      await session.step(137, "When user picks \"Reset View\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Reset View", el("line chart viewer")));
      await session.step(138, "Then \"d4-linechart-reset-view\" event should have fired on line chart viewer", () => eventFired(page, "d4-linechart-reset-view", el("line chart viewer")));
      await session.step(139, "And the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("X Min and X Max pin the window, and Reset View returns to them rather than to the column", async () => {
      await session.step(143, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMin","634800"],["xMax","634850"]]));
      await session.step(146, "Then the \"x axis min\" reading of line chart viewer should be 634800", () => readingIs(page, "x axis min", el("line chart viewer"), 634800));
      await session.step(147, "And the \"x axis max\" reading of line chart viewer should be 634850", () => readingIs(page, "x axis max", el("line chart viewer"), 634850));
      await session.step(148, "And the \"x axis span\" reading of line chart viewer should be 50", () => readingIs(page, "x axis span", el("line chart viewer"), 50));
      await session.step(149, "And the \"markers drawn\" reading of line chart viewer should be lower than before", () => readingLower(page, "markers drawn", el("line chart viewer")));
      await session.step(150, "When user scrolls the mouse wheel up over the \"plot\" area of line chart viewer", () => wheelOverArea(page, "up", "plot", el("line chart viewer")));
      await session.step(151, "Then the \"x axis span\" reading of line chart viewer should be lower than before", () => readingLower(page, "x axis span", el("line chart viewer")));
      await session.step(152, "When user picks \"Reset View\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Reset View", el("line chart viewer")));
      await session.step(153, "Then \"xMin\" property of line chart viewer should be \"634800\"", () => propertyShouldBe(page, "xMin", el("line chart viewer"), "634800"));
      await session.step(154, "And \"xMax\" property of line chart viewer should be \"634850\"", () => propertyShouldBe(page, "xMax", el("line chart viewer"), "634850"));
      await session.step(155, "And the \"x axis span\" reading of line chart viewer should be 50", () => readingIs(page, "x axis span", el("line chart viewer"), 50));
      await session.step(156, "When user sets properties of line chart viewer:", () => setProperties(page, el("line chart viewer"), [["xMin",""],["xMax",""]]));
      await session.step(159, "Then the \"x axis span\" reading of line chart viewer should be between 106 and 107", () => readingBetween(page, "x axis span", el("line chart viewer"), 106, 107));
      await session.step(160, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
