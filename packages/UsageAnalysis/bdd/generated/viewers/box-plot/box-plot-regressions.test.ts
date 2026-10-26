/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, core#307, public#3469, public#2764, public#3066, public#1379]
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
import {addCategoricalFilter, clearSelection, colorConditional, colorLinear, colorOff, filterPasses, noneOfSelected, onlyOfSelected, selectedRowCount, someOfSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, doubleClickArea, dragSelectionOverArea, hasArea, hasNoArea, lessHighlight, lessInk, moreHighlight, moreInk, narrowerRange, noErrors, propertyShouldBe, rememberRange, rememberedRange, setProperty, widerRange} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {zoomValueAxis} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot regressions from closed issues", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-regressions.feature", import.meta.url);
  test("Toggling Show Empty Categories keeps the zoom (core#307)", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:core#307"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(19, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","Average Mass"],["Category 1","Series"]]), [["Value","Average Mass"],["Category 1","Series"]]);
    await session.step(22, "Then \"Show Empty Categories\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Empty Categories", el("box plot viewer"), "true"));
    await session.step(23, "When user zooms into the value axis of box plot viewer", () => zoomValueAxis(page, el("box plot viewer")));
    await session.step(24, "Then box plot viewer should show a narrower value range than before", () => narrowerRange(page, el("box plot viewer")));
    await session.step(25, "When user remembers the value range of box plot viewer", () => rememberRange(page, el("box plot viewer")));
    await session.step(26, "And user sets \"Show Empty Categories\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "false"));
    await session.step(27, "Then box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
    await session.step(28, "When user sets \"Show Empty Categories\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Empty Categories", el("box plot viewer"), "true"));
    await session.step(29, "Then box plot viewer should show the remembered value range", () => rememberedRange(page, el("box plot viewer")));
    await session.step(30, "When user double-clicks on the \"view\" area of box plot viewer", () => doubleClickArea(page, "view", el("box plot viewer")));
    await session.step(31, "Then box plot viewer should show a wider value range than before", () => widerRange(page, el("box plot viewer")));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
  test("An explicit Marker Color Column survives a Category 1 change (public#3469)", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:public#3469"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(37, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(38, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]), [["Value","AGE"],["Category 1","SEX"]]);
    await session.step(41, "And user sets \"Marker Color Column\" property of box plot viewer to \"HEIGHT\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "HEIGHT"));
    await session.step(42, "Then box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
    await session.step(43, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
    await session.step(44, "Then box plot viewer should have a \"category Caucasian\" area", () => hasArea(page, el("box plot viewer"), "category Caucasian"));
    await session.step(45, "And \"Marker Color Column\" property of box plot viewer should be \"HEIGHT\"", () => propertyShouldBe(page, "Marker Color Column", el("box plot viewer"), "HEIGHT"));
    await session.step(46, "And box plot viewer should have a \"color scale\" area", () => hasArea(page, el("box plot viewer"), "color scale"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Shift-drag picks its category when the categories to its left are filtered out (public#2764)", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:public#2764"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(52, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(53, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(56, "And user adds a categorical filter on \"RACE\" keeping \"Caucasian, Other\"", () => addCategoricalFilter(page, "RACE", "Caucasian, Other"));
    await session.step(57, "Then 958 rows should pass the filter", () => filterPasses(page, 958));
    await session.step(58, "And box plot viewer should not have a \"category Asian\" area", () => hasNoArea(page, el("box plot viewer"), "category Asian"));
    await session.step(59, "When user drags a selection box over the \"Other values\" area of box plot viewer", () => dragSelectionOverArea(page, "Other values", el("box plot viewer")));
    await session.step(60, "Then some rows where \"RACE\" is \"Other\" should be selected", () => someOfSelected(page, "RACE", "Other"));
    await session.step(61, "And no rows where \"RACE\" is \"Caucasian\" should be selected", () => noneOfSelected(page, "RACE", "Caucasian"));
    await session.step(62, "And 51 rows should be selected", () => selectedRowCount(page, 51));
    await session.step(63, "When user clears the row selection", () => clearSelection(page));
    await session.step(64, "And user drags a selection box over the \"Caucasian values\" area of box plot viewer", () => dragSelectionOverArea(page, "Caucasian values", el("box plot viewer")));
    await session.step(65, "Then some rows where \"RACE\" is \"Caucasian\" should be selected", () => someOfSelected(page, "RACE", "Caucasian"));
    await session.step(66, "And no rows where \"RACE\" is \"Other\" should be selected", () => noneOfSelected(page, "RACE", "Other"));
    await session.step(67, "And 768 rows should be selected", () => selectedRowCount(page, 768));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
  });
  test("The selection is highlighted under a linear and under a conditional colouring (public#3066)", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:public#3066"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(73, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(74, "When user colors \"HEIGHT\" column linearly from \"#00FF00\" to \"#0000FF\"", () => colorLinear(page, "HEIGHT", "#00FF00", "#0000FF"));
    await session.step(75, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"],["Marker Color Column","HEIGHT"]]), [["Value","AGE"],["Category 1","RACE"],["Marker Color Column","HEIGHT"]]);
    await session.step(79, "And user clicks on the \"category Caucasian\" area of box plot viewer", () => clickArea(page, "category Caucasian", el("box plot viewer")));
    await session.step(80, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
    await session.step(81, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
    await session.step(82, "When user clears the row selection", () => clearSelection(page));
    await session.step(83, "Then box plot viewer should show less selection highlight than before", () => lessHighlight(page, el("box plot viewer")));
    await session.step(84, "When user colors \"HEIGHT\" column conditionally:", () => colorConditional(page, "HEIGHT", [["> 170","#FF0000"]]), [["> 170","#FF0000"]]);
    await session.step(86, "And user clicks on the \"category Caucasian\" area of box plot viewer", () => clickArea(page, "category Caucasian", el("box plot viewer")));
    await session.step(87, "Then only rows where \"RACE\" is \"Caucasian\" should be selected", () => onlyOfSelected(page, "RACE", "Caucasian"));
    await session.step(88, "And box plot viewer should show more selection highlight than before", () => moreHighlight(page, el("box plot viewer")));
    await session.step(89, "When user clears the row selection", () => clearSelection(page));
    await session.step(90, "And user removes the coloring of \"HEIGHT\" column", () => colorOff(page, "HEIGHT"));
    await session.step(91, "Then no errors should have been logged", () => noErrors(page));
  });
  test("Show Outside Values off takes the outlier markers away (public#1379)", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:public#1379"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(96, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(97, "When user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","WEIGHT"],["Category 1","RACE"]]), [["Value","WEIGHT"],["Category 1","RACE"]]);
    await session.step(100, "And user sets \"Show Inside Values\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Inside Values", el("box plot viewer"), "false"));
    await session.step(101, "And user sets \"Show Outside Values\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Outside Values", el("box plot viewer"), "false"));
    await session.step(102, "Then box plot viewer should have less ink than before", () => lessInk(page, el("box plot viewer")));
    await session.step(103, "And box plot viewer should not have a \"marker\" area", () => hasNoArea(page, el("box plot viewer"), "marker"));
    await session.step(104, "When user sets \"Show Outside Values\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Outside Values", el("box plot viewer"), "true"));
    await session.step(105, "Then box plot viewer should have more ink than before", () => moreInk(page, el("box plot viewer")));
    await session.step(106, "And box plot viewer should have a \"marker\" area", () => hasArea(page, el("box plot viewer"), "marker"));
    await session.step(107, "And no errors should have been logged", () => noErrors(page));
  });
});
