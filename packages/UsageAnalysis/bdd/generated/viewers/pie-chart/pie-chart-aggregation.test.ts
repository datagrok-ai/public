/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart-aggregation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {addCalculated, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hasNoArea, hoverArea, lessInk, moreInk, noErrors, painted, pointerAway, propertyShouldBe, readingBetween, readingHigher, readingIs, readingLower, readingReads, readingsDiffer, readingsEqual, repaintedBy, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart aggregations, validation and the date category map", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart-aggregation.feature", import.meta.url);
  test("Pie chart aggregations, validation and the date category map", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]));
    await session.step(17, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(18, "And \"Segment Angle Column\" property of pie chart viewer should be \"AGE\"", () => propertyShouldBe(page, "Segment Angle Column", el("pie chart viewer"), "AGE"));
    await session.step(19, "And \"Segment Angle Aggr Type\" property of pie chart viewer should be \"count\"", () => propertyShouldBe(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
    await session.step(20, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
    await run.scenario("The angle aggregation walked over the list leaves a valid disc every time", async () => {
      await session.step(23, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"min\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "min"));
      await session.step(24, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
      await session.step(25, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(26, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(27, "And the \"angle value of Asian\" reading of pie chart viewer should be 22", () => readingIs(page, "angle value of Asian", el("pie chart viewer"), 22));
      await session.step(28, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(29, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"max\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "max"));
      await session.step(30, "Then the \"angle value of Asian\" reading of pie chart viewer should be 64", () => readingIs(page, "angle value of Asian", el("pie chart viewer"), 64));
      await session.step(31, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(32, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"med\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "med"));
      await session.step(33, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
      await session.step(34, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(35, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"stdev\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "stdev"));
      await session.step(36, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
      await session.step(37, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(38, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(39, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(40, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"count\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
      await session.step(41, "Then the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("avg and sum of the same column give different shares", async () => {
      await session.step(45, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"avg\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "avg"));
      await session.step(46, "Then the \"share of Caucasian\" reading of pie chart viewer should be between 25.7 and 25.75", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 25.7, 25.75));
      await session.step(47, "And the \"share of Asian\" reading of pie chart viewer should be between 21.4 and 21.5", () => readingBetween(page, "share of Asian", el("pie chart viewer"), 21.4, 21.5));
      await session.step(48, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(49, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"sum\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "sum"));
      await session.step(50, "Then the \"share of Caucasian\" reading of pie chart viewer should be between 89.55 and 89.6", () => readingBetween(page, "share of Caucasian", el("pie chart viewer"), 89.55, 89.6));
      await session.step(51, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 40924", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 40924));
      await session.step(52, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(53, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"count\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
      await session.step(54, "Then the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(55, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Segment Length Column shortens every wedge but the longest", async () => {
      await session.step(59, "Then the \"outer radius of Caucasian\" and \"pie radius\" readings of pie chart viewer should be the same", () => readingsEqual(page, "outer radius of Caucasian", "pie radius", el("pie chart viewer")));
      await session.step(60, "And the \"outer radius of Asian\" and \"pie radius\" readings of pie chart viewer should be the same", () => readingsEqual(page, "outer radius of Asian", "pie radius", el("pie chart viewer")));
      await session.step(61, "When user sets \"Segment Length Column\" property of pie chart viewer to \"WEIGHT\"", () => setProperty(page, "Segment Length Column", el("pie chart viewer"), "WEIGHT"));
      await session.step(62, "Then the \"outer radius of Black\" and \"pie radius\" readings of pie chart viewer should be the same", () => readingsEqual(page, "outer radius of Black", "pie radius", el("pie chart viewer")));
      await session.step(63, "And the \"outer radius of Caucasian\" and \"pie radius\" readings of pie chart viewer should differ", () => readingsDiffer(page, "outer radius of Caucasian", "pie radius", el("pie chart viewer")));
      await session.step(64, "And the \"outer radius of Asian\" reading of pie chart viewer should be lower than before", () => readingLower(page, "outer radius of Asian", el("pie chart viewer")));
      await session.step(65, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(66, "And pie chart viewer should have less ink than before", () => lessInk(page, el("pie chart viewer")));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing the length column gives every wedge the full radius back", async () => {
      await session.step(70, "When user sets \"Segment Length Column\" property of pie chart viewer to \"\"", () => setProperty(page, "Segment Length Column", el("pie chart viewer"), ""));
      await session.step(71, "Then the \"outer radius of Caucasian\" and \"pie radius\" readings of pie chart viewer should be the same", () => readingsEqual(page, "outer radius of Caucasian", "pie radius", el("pie chart viewer")));
      await session.step(72, "And the \"outer radius of Asian\" reading of pie chart viewer should be higher than before", () => readingHigher(page, "outer radius of Asian", el("pie chart viewer")));
      await session.step(73, "And pie chart viewer should have more ink than before", () => moreInk(page, el("pie chart viewer")));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A negative aggregation is refused and says why", async () => {
      await session.step(77, "When user adds a calculated column \"NEG_PROBE\" with formula \"${AGE} - 50\"", () => addCalculated(page, "NEG_PROBE", "${AGE} - 50"));
      await session.step(78, "And user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Segment Angle Column","NEG_PROBE"],["Segment Angle Aggr Type","min"]]));
      await session.step(81, "Then the \"error\" reading of pie chart viewer should be \"min(NEG_PROBE) contains negative values\"", () => readingReads(page, "error", el("pie chart viewer"), "min(NEG_PROBE) contains negative values"));
      await session.step(82, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
      await session.step(83, "And pie chart viewer should not have a \"pie\" area", () => hasNoArea(page, el("pie chart viewer"), "pie"));
      await session.step(84, "And pie chart viewer should not have a \"slice Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "slice Caucasian"));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An all-zero aggregation is refused too, and clearing it draws the disc again", async () => {
      await session.step(88, "When user adds a calculated column \"ZERO_PROBE\" with formula \"0\"", () => addCalculated(page, "ZERO_PROBE", "0"));
      await session.step(89, "And user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Segment Angle Column","ZERO_PROBE"],["Segment Angle Aggr Type","sum"]]));
      await session.step(92, "Then the \"error\" reading of pie chart viewer should be \"sum(ZERO_PROBE): all values are 0\"", () => readingReads(page, "error", el("pie chart viewer"), "sum(ZERO_PROBE): all values are 0"));
      await session.step(93, "And the \"slices\" reading of pie chart viewer should be 0", () => readingIs(page, "slices", el("pie chart viewer"), 0));
      await session.step(94, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Segment Angle Column","AGE"],["Segment Angle Aggr Type","count"]]));
      await session.step(97, "Then the \"error\" reading of pie chart viewer should be \"\"", () => readingReads(page, "error", el("pie chart viewer"), ""));
      await session.step(98, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(99, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(100, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(101, "When user removes \"NEG_PROBE\" column", () => removeColumn(page, "NEG_PROBE"));
      await session.step(102, "And user removes \"ZERO_PROBE\" column", () => removeColumn(page, "ZERO_PROBE"));
      await session.step(103, "Then the table should not have a column \"NEG_PROBE\"", () => hasNoColumn(page, "NEG_PROBE"));
      await session.step(104, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The tooltip names the wedge and follows the configured aggregation", async () => {
      await session.step(107, "When user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(108, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
      await session.step(109, "And tooltip should contain text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(110, "And tooltip should contain text \"896 rows\"", () => shouldContainText(page, el("tooltip"), "896 rows"));
      await session.step(111, "When user moves the pointer away from pie chart viewer", () => pointerAway(page, el("pie chart viewer")));
      await session.step(112, "And user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"avg\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "avg"));
      await session.step(113, "And user hovers over the \"slice Caucasian\" area of pie chart viewer", () => hoverArea(page, "slice Caucasian", el("pie chart viewer")));
      await session.step(114, "Then tooltip should contain text \"avg(AGE): 45.67\"", () => shouldContainText(page, el("tooltip"), "avg(AGE): 45.67"));
      await session.step(115, "When user moves the pointer away from pie chart viewer", () => pointerAway(page, el("pie chart viewer")));
      await session.step(116, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
      await session.step(117, "When user sets \"Segment Angle Aggr Type\" property of pie chart viewer to \"count\"", () => setProperty(page, "Segment Angle Aggr Type", el("pie chart viewer"), "count"));
      await session.step(118, "Then the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(119, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Category Map turns one datetime column into years, months and quarters", async () => {
      await session.step(122, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Category","STARTED"],["Category Map","year"]]));
      await session.step(125, "Then the \"slices\" reading of pie chart viewer should be 3", () => readingIs(page, "slices", el("pie chart viewer"), 3));
      await session.step(126, "And pie chart viewer should have a \"slice 1989\" area", () => hasArea(page, el("pie chart viewer"), "slice 1989"));
      await session.step(127, "And pie chart viewer should have a \"slice 1991\" area", () => hasArea(page, el("pie chart viewer"), "slice 1991"));
      await session.step(128, "And the \"angle value of 1989\" reading of pie chart viewer should be 43", () => readingIs(page, "angle value of 1989", el("pie chart viewer"), 43));
      await session.step(129, "When user sets \"Category Map\" property of pie chart viewer to \"month\"", () => setProperty(page, "Category Map", el("pie chart viewer"), "month"));
      await session.step(130, "Then the \"slices\" reading of pie chart viewer should be 12", () => readingIs(page, "slices", el("pie chart viewer"), 12));
      await session.step(131, "And pie chart viewer should have a \"slice January\" area", () => hasArea(page, el("pie chart viewer"), "slice January"));
      await session.step(132, "And pie chart viewer should not have a \"slice 1989\" area", () => hasNoArea(page, el("pie chart viewer"), "slice 1989"));
      await session.step(133, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(134, "When user sets \"Category Map\" property of pie chart viewer to \"quarter\"", () => setProperty(page, "Category Map", el("pie chart viewer"), "quarter"));
      await session.step(135, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(136, "And pie chart viewer should have a \"slice Q1\" area", () => hasArea(page, el("pie chart viewer"), "slice Q1"));
      await session.step(137, "And pie chart viewer should have a \"slice Q4\" area", () => hasArea(page, el("pie chart viewer"), "slice Q4"));
      await session.step(138, "And pie chart viewer should not have a \"slice January\" area", () => hasNoArea(page, el("pie chart viewer"), "slice January"));
      await session.step(139, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(140, "When user sets \"Category\" property of pie chart viewer to \"RACE\"", () => setProperty(page, "Category", el("pie chart viewer"), "RACE"));
      await session.step(141, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(142, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(143, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
