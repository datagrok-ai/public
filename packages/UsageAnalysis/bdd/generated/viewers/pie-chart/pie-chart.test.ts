/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pie-chart/pie-chart.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pie-chart]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {slicesByCategory, slicesByShareAscending, slicesByShareDescending, slicesOffCentre} from '../../../bindings/pie-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColors, areaNarrower, areaRepainted, areaWider, boundTable, hasArea, hasNoArea, lessInk, moreInk, noErrors, painted, propertyShouldBe, readingHigher, readingIs, readingLower, readingSame, repainted, repaintedBy, resizeTo, restoreSize, setProperties, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Pie chart property surface", () => {
  const session = feature(test, "features/viewers/pie-chart/pie-chart.feature", import.meta.url);
  test("Pie chart property surface", {tag: ["@journey", "@viewers", "@realizes:viewers.pie-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 14, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","RACE"]]));
    await session.step(22, "Then the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
    await session.step(23, "And pie chart viewer should show 1000 rows", () => showsRows(page, el("pie chart viewer"), 1000));
    await run.scenario("A wedge per category, sized by its share of the rows", async () => {
      await session.step(26, "Then pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(27, "And pie chart viewer should have a \"slice Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "slice Caucasian"));
      await session.step(28, "And pie chart viewer should have a \"slice Asian\" area", () => hasArea(page, el("pie chart viewer"), "slice Asian"));
      await session.step(29, "And pie chart viewer should have a \"slice Black\" area", () => hasArea(page, el("pie chart viewer"), "slice Black"));
      await session.step(30, "And pie chart viewer should have a \"slice Other\" area", () => hasArea(page, el("pie chart viewer"), "slice Other"));
      await session.step(31, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(32, "And the \"share of Asian\" reading of pie chart viewer should be 1.5", () => readingIs(page, "share of Asian", el("pie chart viewer"), 1.5));
      await session.step(33, "And the \"angle value of Caucasian\" reading of pie chart viewer should be 896", () => readingIs(page, "angle value of Caucasian", el("pie chart viewer"), 896));
      await session.step(34, "And the \"angle value of Asian\" reading of pie chart viewer should be 15", () => readingIs(page, "angle value of Asian", el("pie chart viewer"), 15));
      await session.step(35, "And the \"pie\" area of pie chart viewer should be painted in at least 4 colors", () => areaColors(page, "pie", el("pie chart viewer"), 4));
      await session.step(36, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sorting by value puts the smallest wedge first and the order reverses", async () => {
      await session.step(39, "Then \"Pie Sort Type\" property of pie chart viewer should be \"by value\"", () => propertyShouldBe(page, "Pie Sort Type", el("pie chart viewer"), "by value"));
      await session.step(40, "And the slices of pie chart viewer should be ordered by share ascending", () => slicesByShareAscending(page, el("pie chart viewer")));
      await session.step(41, "And the \"start angle of Asian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 0));
      await session.step(42, "When user sets \"Pie Sort Order\" property of pie chart viewer to \"desc\"", () => setProperty(page, "Pie Sort Order", el("pie chart viewer"), "desc"));
      await session.step(43, "Then the slices of pie chart viewer should be ordered by share descending", () => slicesByShareDescending(page, el("pie chart viewer")));
      await session.step(44, "And the \"start angle of Caucasian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Caucasian", el("pie chart viewer"), 0));
      await session.step(45, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(46, "When user sets \"Pie Sort Order\" property of pie chart viewer to \"asc\"", () => setProperty(page, "Pie Sort Order", el("pie chart viewer"), "asc"));
      await session.step(47, "Then the slices of pie chart viewer should be ordered by share ascending", () => slicesByShareAscending(page, el("pie chart viewer")));
      await session.step(48, "And the \"start angle of Asian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 0));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Sorting by category draws the wedges alphabetically", async () => {
      await session.step(52, "When user sets \"Pie Sort Type\" property of pie chart viewer to \"by category\"", () => setProperty(page, "Pie Sort Type", el("pie chart viewer"), "by category"));
      await session.step(53, "Then the slices of pie chart viewer should be ordered by category", () => slicesByCategory(page, el("pie chart viewer")));
      await session.step(54, "And the \"start angle of Asian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 0));
      await session.step(55, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(56, "When user sets \"Pie Sort Type\" property of pie chart viewer to \"by value\"", () => setProperty(page, "Pie Sort Type", el("pie chart viewer"), "by value"));
      await session.step(57, "Then the slices of pie chart viewer should be ordered by share ascending", () => slicesByShareAscending(page, el("pie chart viewer")));
      await session.step(58, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Start Angle rotates the whole disc", async () => {
      await session.step(62, "Then the \"start angle of Asian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 0));
      await session.step(63, "When user sets \"Start Angle\" property of pie chart viewer to \"90\"", () => setProperty(page, "Start Angle", el("pie chart viewer"), "90"));
      await session.step(64, "Then the \"start angle of Asian\" reading of pie chart viewer should be 90", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 90));
      await session.step(65, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(66, "When user sets \"Start Angle\" property of pie chart viewer to \"180\"", () => setProperty(page, "Start Angle", el("pie chart viewer"), "180"));
      await session.step(67, "Then the \"start angle of Asian\" reading of pie chart viewer should be 180", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 180));
      await session.step(68, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(69, "When user sets \"Start Angle\" property of pie chart viewer to \"0\"", () => setProperty(page, "Start Angle", el("pie chart viewer"), "0"));
      await session.step(70, "Then the \"start angle of Asian\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of Asian", el("pie chart viewer"), 0));
      await session.step(71, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Max Radius caps the disc and gives it back", async () => {
      await session.step(75, "Then the \"pie radius\" reading of pie chart viewer should be 150", () => readingIs(page, "pie radius", el("pie chart viewer"), 150));
      await session.step(76, "When user sets \"Max Radius\" property of pie chart viewer to \"100\"", () => setProperty(page, "Max Radius", el("pie chart viewer"), "100"));
      await session.step(77, "Then the \"pie radius\" reading of pie chart viewer should be 100", () => readingIs(page, "pie radius", el("pie chart viewer"), 100));
      await session.step(78, "And the \"outer radius of Caucasian\" reading of pie chart viewer should be 100", () => readingIs(page, "outer radius of Caucasian", el("pie chart viewer"), 100));
      await session.step(79, "And the \"pie\" area of pie chart viewer should be narrower than before", () => areaNarrower(page, "pie", el("pie chart viewer")));
      await session.step(80, "And pie chart viewer should have less ink than before", () => lessInk(page, el("pie chart viewer")));
      await session.step(81, "When user sets \"Max Radius\" property of pie chart viewer to \"150\"", () => setProperty(page, "Max Radius", el("pie chart viewer"), "150"));
      await session.step(82, "Then the \"pie radius\" reading of pie chart viewer should be 150", () => readingIs(page, "pie radius", el("pie chart viewer"), 150));
      await session.step(83, "And the \"pie\" area of pie chart viewer should be wider than before", () => areaWider(page, "pie", el("pie chart viewer")));
      await session.step(84, "And pie chart viewer should have more ink than before", () => moreInk(page, el("pie chart viewer")));
      await session.step(85, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Shift explodes the wedges out of the centre", async () => {
      await session.step(88, "Then the slices of pie chart viewer should sit 0 pixels off the centre", () => slicesOffCentre(page, el("pie chart viewer"), 0));
      await session.step(89, "When user sets \"Shift\" property of pie chart viewer to \"20\"", () => setProperty(page, "Shift", el("pie chart viewer"), "20"));
      await session.step(90, "Then the slices of pie chart viewer should sit 20 pixels off the centre", () => slicesOffCentre(page, el("pie chart viewer"), 20));
      await session.step(91, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(92, "When user sets \"Shift\" property of pie chart viewer to \"0\"", () => setProperty(page, "Shift", el("pie chart viewer"), "0"));
      await session.step(93, "Then the slices of pie chart viewer should sit 0 pixels off the centre", () => slicesOffCentre(page, el("pie chart viewer"), 0));
      await session.step(94, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Label Position moves the labels out of the wedges and back in", async () => {
      await session.step(98, "Then pie chart viewer should have a \"label of Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "label of Caucasian"));
      await session.step(99, "And pie chart viewer should have a \"label of Asian\" area", () => hasArea(page, el("pie chart viewer"), "label of Asian"));
      await session.step(100, "And the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(101, "When user sets \"Label Position\" property of pie chart viewer to \"Inside\"", () => setProperty(page, "Label Position", el("pie chart viewer"), "Inside"));
      await session.step(102, "Then pie chart viewer should not have a \"label of Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "label of Caucasian"));
      await session.step(103, "And the \"labels shown\" reading of pie chart viewer should be 3", () => readingIs(page, "labels shown", el("pie chart viewer"), 3));
      await session.step(104, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(105, "When user sets \"Label Position\" property of pie chart viewer to \"Outside\"", () => setProperty(page, "Label Position", el("pie chart viewer"), "Outside"));
      await session.step(106, "Then pie chart viewer should have a \"label of Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "label of Caucasian"));
      await session.step(107, "And pie chart viewer should have a \"label of Asian\" area", () => hasArea(page, el("pie chart viewer"), "label of Asian"));
      await session.step(108, "And the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(109, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(110, "When user sets \"Label Position\" property of pie chart viewer to \"Auto\"", () => setProperty(page, "Label Position", el("pie chart viewer"), "Auto"));
      await session.step(111, "Then pie chart viewer should have a \"label of Caucasian\" area", () => hasArea(page, el("pie chart viewer"), "label of Caucasian"));
      await session.step(112, "And the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With nothing left to say a wedge gets no label at all", async () => {
      await session.step(116, "Then the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(117, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Show Label","false"],["Show Percentage","false"]]));
      await session.step(120, "Then the \"labels shown\" reading of pie chart viewer should be 0", () => readingIs(page, "labels shown", el("pie chart viewer"), 0));
      await session.step(121, "And pie chart viewer should not have a \"label of Caucasian\" area", () => hasNoArea(page, el("pie chart viewer"), "label of Caucasian"));
      await session.step(122, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(123, "When user sets \"Show Value\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Value", el("pie chart viewer"), "true"));
      await session.step(124, "Then the \"labels shown\" reading of pie chart viewer should be higher than before", () => readingHigher(page, "labels shown", el("pie chart viewer")));
      await session.step(125, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(126, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Show Label","true"],["Show Percentage","true"],["Show Value","false"]]));
      await session.step(130, "Then the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(131, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Outline Line Width thickens the border between the wedges", async () => {
      await session.step(135, "When user sets \"Outline Line Width\" property of pie chart viewer to \"5\"", () => setProperty(page, "Outline Line Width", el("pie chart viewer"), "5"));
      await session.step(136, "Then \"Outline Line Width\" property of pie chart viewer should be \"5\"", () => propertyShouldBe(page, "Outline Line Width", el("pie chart viewer"), "5"));
      await session.step(137, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(138, "When user sets \"Outline Line Width\" property of pie chart viewer to \"1\"", () => setProperty(page, "Outline Line Width", el("pie chart viewer"), "1"));
      await session.step(139, "Then \"Outline Line Width\" property of pie chart viewer should be \"1\"", () => propertyShouldBe(page, "Outline Line Width", el("pie chart viewer"), "1"));
      await session.step(140, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(141, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Column Selector takes the on-chart selector off the canvas", async () => {
      await session.step(144, "Then pie chart viewer should have a \"column selector\" area", () => hasArea(page, el("pie chart viewer"), "column selector"));
      await session.step(145, "When user sets \"Show Column Selector\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Column Selector", el("pie chart viewer"), "false"));
      await session.step(146, "Then pie chart viewer should not have a \"column selector\" area", () => hasNoArea(page, el("pie chart viewer"), "column selector"));
      await session.step(147, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(148, "When user sets \"Show Column Selector\" property of pie chart viewer to \"true\"", () => setProperty(page, "Show Column Selector", el("pie chart viewer"), "true"));
      await session.step(149, "Then pie chart viewer should have a \"column selector\" area", () => hasArea(page, el("pie chart viewer"), "column selector"));
      await session.step(150, "And pie chart viewer should have repainted by at least 300 pixels", () => repaintedBy(page, el("pie chart viewer"), 300));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Auto Layout off hands the margins to the chart box", async () => {
      await session.step(154, "Then the \"pie radius\" reading of pie chart viewer should be 150", () => readingIs(page, "pie radius", el("pie chart viewer"), 150));
      await session.step(155, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Auto Layout","false"],["Margin Left","100"],["Margin Top","100"]]));
      await session.step(159, "Then the \"pie radius\" reading of pie chart viewer should be lower than before", () => readingLower(page, "pie radius", el("pie chart viewer")));
      await session.step(160, "And the \"view\" area of pie chart viewer should be narrower than before", () => areaNarrower(page, "view", el("pie chart viewer")));
      await session.step(161, "And pie chart viewer should have repainted by at least 500 pixels", () => repaintedBy(page, el("pie chart viewer"), 500));
      await session.step(162, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Auto Layout","true"],["Margin Left","10"],["Margin Top","10"]]));
      await session.step(166, "Then the \"pie radius\" reading of pie chart viewer should be 150", () => readingIs(page, "pie radius", el("pie chart viewer"), 150));
      await session.step(167, "And the \"view\" area of pie chart viewer should be wider than before", () => areaWider(page, "view", el("pie chart viewer")));
      await session.step(168, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A viewer too small for its labels drops them and the selector", async () => {
      await session.step(171, "Then the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(172, "And pie chart viewer should have a \"column selector\" area", () => hasArea(page, el("pie chart viewer"), "column selector"));
      await session.step(173, "When user resizes pie chart viewer to 150 by 150", () => resizeTo(page, el("pie chart viewer"), 150, 150));
      await session.step(174, "Then the \"labels shown\" reading of pie chart viewer should be 0", () => readingIs(page, "labels shown", el("pie chart viewer"), 0));
      await session.step(175, "And pie chart viewer should not have a \"column selector\" area", () => hasNoArea(page, el("pie chart viewer"), "column selector"));
      await session.step(176, "And the \"pie radius\" reading of pie chart viewer should be lower than before", () => readingLower(page, "pie radius", el("pie chart viewer")));
      await session.step(177, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(178, "When user restores the size of pie chart viewer", () => restoreSize(page, el("pie chart viewer")));
      await session.step(179, "Then the \"labels shown\" reading of pie chart viewer should be 4", () => readingIs(page, "labels shown", el("pie chart viewer"), 4));
      await session.step(180, "And pie chart viewer should have a \"column selector\" area", () => hasArea(page, el("pie chart viewer"), "column selector"));
      await session.step(181, "And the \"pie radius\" reading of pie chart viewer should be 150", () => readingIs(page, "pie radius", el("pie chart viewer"), 150));
      await session.step(182, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Table property rebinds the chart to another table", async () => {
      await session.step(185, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(186, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(187, "When user sets \"Table\" property of pie chart viewer to \"spgi-100\"", () => setProperty(page, "Table", el("pie chart viewer"), "spgi-100"));
      await session.step(188, "Then pie chart viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("pie chart viewer"), "spgi-100"));
      await session.step(189, "When user sets \"Category\" property of pie chart viewer to \"Primary Series Name\"", () => setProperty(page, "Category", el("pie chart viewer"), "Primary Series Name"));
      await session.step(190, "Then the \"slices\" reading of pie chart viewer should be 5", () => readingIs(page, "slices", el("pie chart viewer"), 5));
      await session.step(191, "And pie chart viewer should show 100 rows", () => showsRows(page, el("pie chart viewer"), 100));
      await session.step(192, "And pie chart viewer should be painted", () => painted(page, el("pie chart viewer")));
      await session.step(193, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Table","demog-1000"],["Category","RACE"]]));
      await session.step(196, "Then pie chart viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("pie chart viewer"), "demog-1000"));
      await session.step(197, "And the \"slices\" reading of pie chart viewer should be 4", () => readingIs(page, "slices", el("pie chart viewer"), 4));
      await session.step(198, "And pie chart viewer should show 1000 rows", () => showsRows(page, el("pie chart viewer"), 1000));
      await session.step(199, "And the \"share of Caucasian\" reading of pie chart viewer should be 89.6", () => readingIs(page, "share of Caucasian", el("pie chart viewer"), 89.6));
      await session.step(200, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Donut mode opens a hole and puts the column name in it", async () => {
      await session.step(203, "Then pie chart viewer should not have a \"donut hole\" area", () => hasNoArea(page, el("pie chart viewer"), "donut hole"));
      await session.step(204, "And pie chart viewer should not have a \"centre label\" area", () => hasNoArea(page, el("pie chart viewer"), "centre label"));
      await session.step(205, "When user sets \"Mode\" property of pie chart viewer to \"Donut\"", () => setProperty(page, "Mode", el("pie chart viewer"), "Donut"));
      await session.step(206, "Then pie chart viewer should have a \"donut hole\" area", () => hasArea(page, el("pie chart viewer"), "donut hole"));
      await session.step(207, "And pie chart viewer should have a \"centre label\" area", () => hasArea(page, el("pie chart viewer"), "centre label"));
      await session.step(208, "And the \"pie radius\" reading of pie chart viewer should be the same as before", () => readingSame(page, "pie radius", el("pie chart viewer")));
      await session.step(209, "And pie chart viewer should have less ink than before", () => lessInk(page, el("pie chart viewer")));
      await session.step(210, "And pie chart viewer should have repainted by at least 1000 pixels", () => repaintedBy(page, el("pie chart viewer"), 1000));
      await session.step(211, "When user sets \"Center Label\" property of pie chart viewer to \"Race mix\"", () => setProperty(page, "Center Label", el("pie chart viewer"), "Race mix"));
      await session.step(212, "Then the \"centre label\" area of pie chart viewer should have repainted", () => areaRepainted(page, "centre label", el("pie chart viewer")));
      await session.step(213, "When user sets \"Show Center Label\" property of pie chart viewer to \"false\"", () => setProperty(page, "Show Center Label", el("pie chart viewer"), "false"));
      await session.step(214, "Then pie chart viewer should not have a \"centre label\" area", () => hasNoArea(page, el("pie chart viewer"), "centre label"));
      await session.step(215, "And pie chart viewer should have a \"donut hole\" area", () => hasArea(page, el("pie chart viewer"), "donut hole"));
      await session.step(216, "And pie chart viewer should have repainted", () => repainted(page, el("pie chart viewer")));
      await session.step(217, "When user sets properties of pie chart viewer:", () => setProperties(page, el("pie chart viewer"), [["Show Center Label","true"],["Center Label",""],["Mode","Pie"]]));
      await session.step(221, "Then pie chart viewer should not have a \"donut hole\" area", () => hasNoArea(page, el("pie chart viewer"), "donut hole"));
      await session.step(222, "And pie chart viewer should not have a \"centre label\" area", () => hasNoArea(page, el("pie chart viewer"), "centre label"));
      await session.step(223, "And pie chart viewer should have more ink than before", () => moreInk(page, el("pie chart viewer")));
      await session.step(224, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
