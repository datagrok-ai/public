/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotation-regions/annotation-regions-interaction-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot, viewers.density-plot, viewers.histogram, viewers.bar-chart]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, noneSelected, onlyBetweenSelected, onlyOfAnySelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areasSameSize, clickArea, clickAreaHolding, hasArea, hoverArea, noErrors, pointerAway, readingIs, resizeTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Hovering and clicking regions on every viewer", () => {
  const session = feature(test, "features/annotation-regions/annotation-regions-interaction-viewers.feature", import.meta.url);
  test("Two overlapping bands on the scatter plot", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "Given user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["markerDefaultSize","2"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Tall\",\"formula1\":\"${HEIGHT} = 180\",\"formula2\":\"${HEIGHT} = 200\"},{\"type\":\"formula\",\"header\":\"Heavy\",\"formula1\":\"${WEIGHT} = 110\",\"formula2\":\"${WEIGHT} = 170\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["markerDefaultSize","2"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Tall\",\"formula1\":\"${HEIGHT} = 180\",\"formula2\":\"${HEIGHT} = 200\"},{\"type\":\"formula\",\"header\":\"Heavy\",\"formula1\":\"${WEIGHT} = 110\",\"formula2\":\"${WEIGHT} = 170\"}]"]]);
    await session.step(28, "And user resizes scatter plot viewer to 800 by 500", () => resizeTo(page, el("scatter plot viewer"), 800, 500));
    await session.step(29, "Then the \"regions shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions shown", el("scatter plot viewer"), 2));
    await session.step(30, "And the \"region titles shown\" reading of scatter plot viewer should be 2", () => readingIs(page, "region titles shown", el("scatter plot viewer"), 2));
    await session.step(31, "When user hovers over the \"left edge of region Tall\" area of scatter plot viewer", () => hoverArea(page, "left edge of region Tall", el("scatter plot viewer")));
    await session.step(32, "Then the \"regions hovered\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 1));
    await session.step(33, "And tooltip should contain text \"Tall\"", () => shouldContainText(page, el("tooltip"), "Tall"));
    await session.step(34, "And tooltip should contain text \"159 rows\"", () => shouldContainText(page, el("tooltip"), "159 rows"));
    await session.step(35, "When user hovers over the \"left edge of overlap of region Tall and region Heavy\" area of scatter plot viewer", () => hoverArea(page, "left edge of overlap of region Tall and region Heavy", el("scatter plot viewer")));
    await session.step(36, "Then the \"regions hovered\" reading of scatter plot viewer should be 2", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 2));
    await session.step(37, "And tooltip should contain text \"Tall\"", () => shouldContainText(page, el("tooltip"), "Tall"));
    await session.step(38, "And tooltip should contain text \"Heavy\"", () => shouldContainText(page, el("tooltip"), "Heavy"));
    await session.step(39, "And tooltip should contain text \"21 rows\"", () => shouldContainText(page, el("tooltip"), "21 rows"));
    await session.step(40, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(41, "Then the \"regions hovered\" reading of scatter plot viewer should be 0", () => readingIs(page, "regions hovered", el("scatter plot viewer"), 0));
    await session.step(42, "When user hovers over the \"left edge of region Tall\" area of scatter plot viewer", () => hoverArea(page, "left edge of region Tall", el("scatter plot viewer")));
    await session.step(43, "And user clicks on the \"left edge of region Tall\" area of scatter plot viewer", () => clickArea(page, "left edge of region Tall", el("scatter plot viewer")));
    await session.step(44, "Then only rows where \"HEIGHT\" is between 180 and 200 should be selected", () => onlyBetweenSelected(page, "HEIGHT", 180, 200));
    await session.step(45, "When user clicks on the \"left edge of region Tall\" area of scatter plot viewer holding Control", () => clickAreaHolding(page, "left edge of region Tall", el("scatter plot viewer"), "Control"));
    await session.step(46, "Then no rows should be selected", () => noneSelected(page));
    await session.step(47, "When user hovers over the \"left edge of overlap of region Tall and region Heavy\" area of scatter plot viewer", () => hoverArea(page, "left edge of overlap of region Tall and region Heavy", el("scatter plot viewer")));
    await session.step(48, "And user clicks on the \"left edge of overlap of region Tall and region Heavy\" area of scatter plot viewer", () => clickArea(page, "left edge of overlap of region Tall and region Heavy", el("scatter plot viewer")));
    await session.step(49, "Then 21 rows should be selected", () => selectedRowCount(page, 21));
    await session.step(50, "When user clears the row selection", () => clearSelection(page));
    await session.step(51, "And user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(52, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Two overlapping bands on the density plot", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "Given user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Tall\",\"formula1\":\"${HEIGHT} = 180\",\"formula2\":\"${HEIGHT} = 200\"},{\"type\":\"formula\",\"header\":\"Heavy\",\"formula1\":\"${WEIGHT} = 110\",\"formula2\":\"${WEIGHT} = 170\"}]"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Tall\",\"formula1\":\"${HEIGHT} = 180\",\"formula2\":\"${HEIGHT} = 200\"},{\"type\":\"formula\",\"header\":\"Heavy\",\"formula1\":\"${WEIGHT} = 110\",\"formula2\":\"${WEIGHT} = 170\"}]"]]);
    await session.step(59, "And user resizes density plot viewer to 800 by 500", () => resizeTo(page, el("density plot viewer"), 800, 500));
    await session.step(60, "Then the \"regions shown\" reading of density plot viewer should be 2", () => readingIs(page, "regions shown", el("density plot viewer"), 2));
    await session.step(61, "When user hovers over the \"left edge of region Tall\" area of density plot viewer", () => hoverArea(page, "left edge of region Tall", el("density plot viewer")));
    await session.step(62, "Then the \"regions hovered\" reading of density plot viewer should be 1", () => readingIs(page, "regions hovered", el("density plot viewer"), 1));
    await session.step(63, "When user hovers over the \"left edge of overlap of region Tall and region Heavy\" area of density plot viewer", () => hoverArea(page, "left edge of overlap of region Tall and region Heavy", el("density plot viewer")));
    await session.step(64, "Then the \"regions hovered\" reading of density plot viewer should be 2", () => readingIs(page, "regions hovered", el("density plot viewer"), 2));
    await session.step(65, "When user clicks on the \"left edge of overlap of region Tall and region Heavy\" area of density plot viewer", () => clickArea(page, "left edge of overlap of region Tall and region Heavy", el("density plot viewer")));
    await session.step(66, "Then 21 rows should be selected", () => selectedRowCount(page, 21));
    await session.step(67, "When user hovers over the \"left edge of region Tall\" area of density plot viewer", () => hoverArea(page, "left edge of region Tall", el("density plot viewer")));
    await session.step(68, "And user clicks on the \"left edge of region Tall\" area of density plot viewer", () => clickArea(page, "left edge of region Tall", el("density plot viewer")));
    await session.step(69, "Then only rows where \"HEIGHT\" is between 180 and 200 should be selected", () => onlyBetweenSelected(page, "HEIGHT", 180, 200));
    await session.step(70, "When user clicks on the \"left edge of region Tall\" area of density plot viewer holding Control", () => clickAreaHolding(page, "left edge of region Tall", el("density plot viewer"), "Control"));
    await session.step(71, "Then no rows should be selected", () => noneSelected(page));
    await session.step(72, "When user moves the pointer away from density plot viewer", () => pointerAway(page, el("density plot viewer")));
    await session.step(73, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A value band on the histogram reacts where no bar is under it", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(76, "Given user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Older\",\"formula1\":\"${AGE} = 60\",\"formula2\":\"${AGE} = 90\"}]"]]), [["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Older\",\"formula1\":\"${AGE} = 60\",\"formula2\":\"${AGE} = 90\"}]"]]);
    await session.step(79, "And user resizes histogram viewer to 800 by 500", () => resizeTo(page, el("histogram viewer"), 800, 500));
    await session.step(80, "Then the \"regions shown\" reading of histogram viewer should be 1", () => readingIs(page, "regions shown", el("histogram viewer"), 1));
    await session.step(81, "When user hovers over the \"top left corner of region Older\" area of histogram viewer", () => hoverArea(page, "top left corner of region Older", el("histogram viewer")));
    await session.step(82, "Then the \"regions hovered\" reading of histogram viewer should be 1", () => readingIs(page, "regions hovered", el("histogram viewer"), 1));
    await session.step(83, "And tooltip should contain text \"Older\"", () => shouldContainText(page, el("tooltip"), "Older"));
    await session.step(84, "When user clicks on the \"top left corner of region Older\" area of histogram viewer", () => clickArea(page, "top left corner of region Older", el("histogram viewer")));
    await session.step(85, "Then only rows where \"AGE\" is between 60 and 90 should be selected", () => onlyBetweenSelected(page, "AGE", 60, 90));
    await session.step(86, "When user clicks on the \"top left corner of region Older\" area of histogram viewer holding Control", () => clickAreaHolding(page, "top left corner of region Older", el("histogram viewer"), "Control"));
    await session.step(87, "Then no rows should be selected", () => noneSelected(page));
    await session.step(88, "When user moves the pointer away from histogram viewer", () => pointerAway(page, el("histogram viewer")));
    await session.step(89, "Then the \"regions hovered\" reading of histogram viewer should be 0", () => readingIs(page, "regions hovered", el("histogram viewer"), 0));
    await session.step(90, "And no errors should have been logged", () => noErrors(page));
  });
  test("A value band on the bar chart reacts between the bars", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(93, "Given user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","vertical"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Mid\",\"formula1\":\"${AGE} = 45.7\",\"formula2\":\"${AGE} = 47.1\"}]"]]), [["splitColumnName","RACE"],["valueColumnName","AGE"],["valueAggrType","avg"],["orientation","vertical"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Mid\",\"formula1\":\"${AGE} = 45.7\",\"formula2\":\"${AGE} = 47.1\"}]"]]);
    await session.step(99, "And user resizes bar chart viewer to 800 by 500", () => resizeTo(page, el("bar chart viewer"), 800, 500));
    await session.step(100, "Then the \"regions shown\" reading of bar chart viewer should be 1", () => readingIs(page, "regions shown", el("bar chart viewer"), 1));
    await session.step(101, "And bar chart viewer should have a \"region Mid\" area", () => hasArea(page, el("bar chart viewer"), "region Mid"));
    await session.step(102, "And the \"region Mid\" and \"view\" areas of bar chart viewer should be the same width", () => areasSameSize(page, "region Mid", "view", el("bar chart viewer"), "width"));
    await session.step(103, "When user hovers over the \"region Mid\" area of bar chart viewer", () => hoverArea(page, "region Mid", el("bar chart viewer")));
    await session.step(104, "Then the \"regions hovered\" reading of bar chart viewer should be 1", () => readingIs(page, "regions hovered", el("bar chart viewer"), 1));
    await session.step(105, "And tooltip should contain text \"Mid\"", () => shouldContainText(page, el("tooltip"), "Mid"));
    await session.step(106, "When user clicks on the \"region Mid\" area of bar chart viewer", () => clickArea(page, "region Mid", el("bar chart viewer")));
    await session.step(107, "Then only rows where \"RACE\" is one of \"Black, Other\" should be selected", () => onlyOfAnySelected(page, "RACE", "Black, Other"));
    await session.step(108, "When user clicks on the \"region Mid\" area of bar chart viewer holding Control", () => clickAreaHolding(page, "region Mid", el("bar chart viewer"), "Control"));
    await session.step(109, "Then no rows should be selected", () => noneSelected(page));
    await session.step(110, "When user moves the pointer away from bar chart viewer", () => pointerAway(page, el("bar chart viewer")));
    await session.step(111, "Then no errors should have been logged", () => noErrors(page));
  });
  test("A value band on the box plot reacts between the boxes", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(114, "Given user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["category1ColumnName","RACE"],["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Middle age\",\"formula1\":\"${AGE} = 36.0\",\"formula2\":\"${AGE} = 56.0\"}]"]]), [["category1ColumnName","RACE"],["valueColumnName","AGE"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Middle age\",\"formula1\":\"${AGE} = 36.0\",\"formula2\":\"${AGE} = 56.0\"}]"]]);
    await session.step(118, "And user resizes box plot viewer to 800 by 500", () => resizeTo(page, el("box plot viewer"), 800, 500));
    await session.step(119, "Then the \"regions shown\" reading of box plot viewer should be 1", () => readingIs(page, "regions shown", el("box plot viewer"), 1));
    await session.step(120, "And box plot viewer should have a \"region Middle age\" area", () => hasArea(page, el("box plot viewer"), "region Middle age"));
    await session.step(121, "And the \"region Middle age\" and \"view\" areas of box plot viewer should be the same width", () => areasSameSize(page, "region Middle age", "view", el("box plot viewer"), "width"));
    await session.step(122, "When user hovers over the \"left edge of region Middle age\" area of box plot viewer", () => hoverArea(page, "left edge of region Middle age", el("box plot viewer")));
    await session.step(123, "Then the \"regions hovered\" reading of box plot viewer should be 1", () => readingIs(page, "regions hovered", el("box plot viewer"), 1));
    await session.step(124, "And tooltip should contain text \"Middle age\"", () => shouldContainText(page, el("tooltip"), "Middle age"));
    await session.step(125, "And tooltip should contain text \"517 rows\"", () => shouldContainText(page, el("tooltip"), "517 rows"));
    await session.step(126, "When user clicks on the \"left edge of region Middle age\" area of box plot viewer", () => clickArea(page, "left edge of region Middle age", el("box plot viewer")));
    await session.step(127, "Then only rows where \"AGE\" is between 36 and 56 should be selected", () => onlyBetweenSelected(page, "AGE", 36, 56));
    await session.step(128, "When user clicks on the \"left edge of region Middle age\" area of box plot viewer holding Control", () => clickAreaHolding(page, "left edge of region Middle age", el("box plot viewer"), "Control"));
    await session.step(129, "Then no rows should be selected", () => noneSelected(page));
    await session.step(130, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
    await session.step(131, "Then the \"regions hovered\" reading of box plot viewer should be 0", () => readingIs(page, "regions hovered", el("box plot viewer"), 0));
    await session.step(132, "And no errors should have been logged", () => noErrors(page));
  });
  test("Two overlapping bands on the line chart", {tag: ["@viewers", "@realizes:viewers.scatter-plot", "@realizes:viewers.density-plot", "@realizes:viewers.histogram", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(135, "Given user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Medium height\",\"formula1\":\"${avg(HEIGHT)} = 165\",\"formula2\":\"${avg(HEIGHT)} = 175\"},{\"type\":\"formula\",\"header\":\"Older\",\"formula1\":\"${AGE} = 60\",\"formula2\":\"${AGE} = 90\"}]"]]), [["xColumnName","AGE"],["yColumnNames","HEIGHT"],["annotationRegions","[{\"type\":\"formula\",\"header\":\"Medium height\",\"formula1\":\"${avg(HEIGHT)} = 165\",\"formula2\":\"${avg(HEIGHT)} = 175\"},{\"type\":\"formula\",\"header\":\"Older\",\"formula1\":\"${AGE} = 60\",\"formula2\":\"${AGE} = 90\"}]"]]);
    await session.step(139, "And user resizes line chart viewer to 800 by 500", () => resizeTo(page, el("line chart viewer"), 800, 500));
    await session.step(140, "Then the \"regions shown\" reading of line chart viewer should be 2", () => readingIs(page, "regions shown", el("line chart viewer"), 2));
    await session.step(141, "When user hovers over the \"left edge of region Medium height\" area of line chart viewer", () => hoverArea(page, "left edge of region Medium height", el("line chart viewer")));
    await session.step(142, "Then the \"regions hovered\" reading of line chart viewer should be 1", () => readingIs(page, "regions hovered", el("line chart viewer"), 1));
    await session.step(143, "And tooltip should contain text \"Medium height\"", () => shouldContainText(page, el("tooltip"), "Medium height"));
    await session.step(144, "And tooltip should contain text \"822 rows\"", () => shouldContainText(page, el("tooltip"), "822 rows"));
    await session.step(145, "When user hovers over the \"left edge of overlap of region Medium height and region Older\" area of line chart viewer", () => hoverArea(page, "left edge of overlap of region Medium height and region Older", el("line chart viewer")));
    await session.step(146, "Then the \"regions hovered\" reading of line chart viewer should be 2", () => readingIs(page, "regions hovered", el("line chart viewer"), 2));
    await session.step(147, "And tooltip should contain text \"Medium height\"", () => shouldContainText(page, el("tooltip"), "Medium height"));
    await session.step(148, "And tooltip should contain text \"Older\"", () => shouldContainText(page, el("tooltip"), "Older"));
    await session.step(149, "And tooltip should contain text \"102 rows\"", () => shouldContainText(page, el("tooltip"), "102 rows"));
    await session.step(150, "When user clicks on the \"left edge of region Medium height\" area of line chart viewer", () => clickArea(page, "left edge of region Medium height", el("line chart viewer")));
    await session.step(151, "Then 822 rows should be selected", () => selectedRowCount(page, 822));
    await session.step(152, "When user clicks on the \"left edge of region Medium height\" area of line chart viewer holding Control", () => clickAreaHolding(page, "left edge of region Medium height", el("line chart viewer"), "Control"));
    await session.step(153, "Then no rows should be selected", () => noneSelected(page));
    await session.step(154, "When user clicks on the \"left edge of overlap of region Medium height and region Older\" area of line chart viewer", () => clickArea(page, "left edge of overlap of region Medium height and region Older", el("line chart viewer")));
    await session.step(155, "Then 102 rows should be selected", () => selectedRowCount(page, 102));
    await session.step(156, "When user clears the row selection", () => clearSelection(page));
    await session.step(157, "And user moves the pointer away from line chart viewer", () => pointerAway(page, el("line chart viewer")));
    await session.step(158, "Then the \"regions hovered\" reading of line chart viewer should be 0", () => readingIs(page, "regions hovered", el("line chart viewer"), 0));
    await session.step(159, "And no errors should have been logged", () => noErrors(page));
  });
});
