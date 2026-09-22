/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-across-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, categoricalColorIs, colorCategorical, colorCodedCategorically, colorOff, noColorCoding, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, openProject, saveAsProject} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, areaNotColor, dragAreaBy, hasArea, legendDocked, legendItemColor, legendItemsDiffer, legendLists, loadLayout, noErrors, pickColorSwatch, pickFromAreaContextMenu, propertyShouldBe, readingReads, saveLayoutToServer, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("One legend column across seven viewers", () => {
  const session = feature(test, "features/viewers/legend/legend-across-viewers.feature", import.meta.url);
  test("One legend column across seven viewers", {tag: ["@journey", "@viewers", "@realizes:viewers.legend"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 15, page);
    await session.step(41, "Given user is logged in", () => loggedIn(page));
    await session.step(42, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(43, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["colorColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnName","WEIGHT"],["yColumnName","HEIGHT"],["colorColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(49, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["valueColumnName","AGE"],["splitColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["valueColumnName","AGE"],["splitColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(54, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["splitColumnNames","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnName","AGE"],["yColumnNames","WEIGHT"],["splitColumnNames","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(60, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","SEX"],["stackColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["splitColumnName","SEX"],["stackColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(65, "And user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["categoryColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["categoryColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(69, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["xColumnNames","RACE"],["Viewer Type","Scatter plot"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["xColumnNames","RACE"],["Viewer Type","Scatter plot"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(74, "And user sets \"colorColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "colorColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(75, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["markerColorColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["categoryColumnNames","RACE"],["valueColumnName","AGE"],["markerColorColumnName","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(81, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
    await session.step(82, "And no errors should have been logged", () => noErrors(page));
    await run.scenario("The scatter plot legend lists the four races in the colors its canvas paints them [viewer=scatter plot]", async () => {
      await session.step(85, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(86, "And the legend of scatter plot viewer should be docked", () => legendDocked(page, el("scatter plot viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of scatter plot viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("scatter plot viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of scatter plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("scatter plot viewer")));
      await session.step(90, "And the \"view\" area of scatter plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("scatter plot viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The histogram legend lists the four races in the colors its canvas paints them [viewer=histogram]", async () => {
      await session.step(85, "Then the legend of histogram viewer should list 4 items", () => legendLists(page, el("histogram viewer"), 4));
      await session.step(86, "And the legend of histogram viewer should be docked", () => legendDocked(page, el("histogram viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of histogram viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("histogram viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of histogram viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("histogram viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of histogram viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("histogram viewer")));
      await session.step(90, "And the \"view\" area of histogram viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("histogram viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The line chart legend lists the four races in the colors its canvas paints them [viewer=line chart]", async () => {
      await session.step(85, "Then the legend of line chart viewer should list 4 items", () => legendLists(page, el("line chart viewer"), 4));
      await session.step(86, "And the legend of line chart viewer should be docked", () => legendDocked(page, el("line chart viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of line chart viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("line chart viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of line chart viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("line chart viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of line chart viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("line chart viewer")));
      await session.step(90, "And the \"view\" area of line chart viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("line chart viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The bar chart legend lists the four races in the colors its canvas paints them [viewer=bar chart]", async () => {
      await session.step(85, "Then the legend of bar chart viewer should list 4 items", () => legendLists(page, el("bar chart viewer"), 4));
      await session.step(86, "And the legend of bar chart viewer should be docked", () => legendDocked(page, el("bar chart viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of bar chart viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("bar chart viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of bar chart viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("bar chart viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of bar chart viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("bar chart viewer")));
      await session.step(90, "And the \"view\" area of bar chart viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("bar chart viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The pie chart legend lists the four races in the colors its canvas paints them [viewer=pie chart]", async () => {
      await session.step(85, "Then the legend of pie chart viewer should list 4 items", () => legendLists(page, el("pie chart viewer"), 4));
      await session.step(86, "And the legend of pie chart viewer should be docked", () => legendDocked(page, el("pie chart viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of pie chart viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("pie chart viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of pie chart viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("pie chart viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of pie chart viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("pie chart viewer")));
      await session.step(90, "And the \"view\" area of pie chart viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("pie chart viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The trellis plot legend lists the four races in the colors its canvas paints them [viewer=trellis plot]", async () => {
      await session.step(85, "Then the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(86, "And the legend of trellis plot viewer should be docked", () => legendDocked(page, el("trellis plot viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of trellis plot viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("trellis plot viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of trellis plot viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("trellis plot viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of trellis plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("trellis plot viewer")));
      await session.step(90, "And the \"view\" area of trellis plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("trellis plot viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The box plot legend lists the four races in the colors its canvas paints them [viewer=box plot]", async () => {
      await session.step(85, "Then the legend of box plot viewer should list 4 items", () => legendLists(page, el("box plot viewer"), 4));
      await session.step(86, "And the legend of box plot viewer should be docked", () => legendDocked(page, el("box plot viewer")));
      await session.step(87, "And the \"Caucasian\" item in the legend of box plot viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("box plot viewer"), "#2CA02C"));
      await session.step(88, "And the \"Other\" item in the legend of box plot viewer should be colored \"#D62728\"", () => legendItemColor(page, "Other", el("box plot viewer"), "#D62728"));
      await session.step(89, "And the \"Caucasian\" and \"Asian\" items in the legend of box plot viewer should be colored differently", () => legendItemsDiffer(page, "Caucasian", "Asian", el("box plot viewer")));
      await session.step(90, "And the \"view\" area of box plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "view", el("box plot viewer"), "#2CA02C"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Other's red is on every canvas that draws the rows itself", async () => {
      await session.step(104, "Then the \"view\" area of scatter plot viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("scatter plot viewer"), "#D62728"));
      await session.step(105, "And the \"view\" area of histogram viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("histogram viewer"), "#D62728"));
      await session.step(106, "And the \"view\" area of line chart viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("line chart viewer"), "#D62728"));
      await session.step(107, "And the \"view\" area of bar chart viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("bar chart viewer"), "#D62728"));
      await session.step(108, "And the \"view\" area of pie chart viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("pie chart viewer"), "#D62728"));
      await session.step(109, "And the \"view\" area of box plot viewer should contain the color \"#D62728\"", () => areaColor(page, "view", el("box plot viewer"), "#D62728"));
      await session.step(110, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A categorical coloring switched on in the grid reaches every legend and every canvas", async () => {
      await session.step(113, "When user drags the \"x scroll handle\" area of grid by 60 pixels to the right", () => dragAreaBy(page, "x scroll handle", el("grid"), 60, "right"));
      await session.step(114, "Then grid should have a \"header RACE\" area", () => hasArea(page, el("grid"), "header RACE"));
      await session.step(115, "When user picks \"Color Coding > Categorical\" from the context menu of the \"header RACE\" area of grid", () => pickFromAreaContextMenu(page, "Color Coding > Categorical", "header RACE", el("grid")));
      await session.step(116, "Then \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(117, "When user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#FF00FF"],["Other","#FFFF00"]]), [["Caucasian","#FF00FF"],["Other","#FFFF00"]]);
      await session.step(120, "Then the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF00FF"));
      await session.step(121, "And the \"Caucasian\" item in the legend of histogram viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("histogram viewer"), "#FF00FF"));
      await session.step(122, "And the \"Caucasian\" item in the legend of line chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("line chart viewer"), "#FF00FF"));
      await session.step(123, "And the \"Caucasian\" item in the legend of bar chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("bar chart viewer"), "#FF00FF"));
      await session.step(124, "And the \"Caucasian\" item in the legend of pie chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("pie chart viewer"), "#FF00FF"));
      await session.step(125, "And the \"Caucasian\" item in the legend of trellis plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("trellis plot viewer"), "#FF00FF"));
      await session.step(126, "And the \"Caucasian\" item in the legend of box plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("box plot viewer"), "#FF00FF"));
      await session.step(127, "And the \"Other\" item in the legend of scatter plot viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("scatter plot viewer"), "#FFFF00"));
      await session.step(128, "And the \"Other\" item in the legend of histogram viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("histogram viewer"), "#FFFF00"));
      await session.step(129, "And the \"Other\" item in the legend of line chart viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("line chart viewer"), "#FFFF00"));
      await session.step(130, "And the \"Other\" item in the legend of bar chart viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("bar chart viewer"), "#FFFF00"));
      await session.step(131, "And the \"Other\" item in the legend of pie chart viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("pie chart viewer"), "#FFFF00"));
      await session.step(132, "And the \"Other\" item in the legend of trellis plot viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("trellis plot viewer"), "#FFFF00"));
      await session.step(133, "And the \"Other\" item in the legend of box plot viewer should be colored \"#FFFF00\"", () => legendItemColor(page, "Other", el("box plot viewer"), "#FFFF00"));
      await session.step(134, "And the \"view\" area of scatter plot viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("scatter plot viewer"), "#FF00FF"));
      await session.step(135, "And the \"view\" area of histogram viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("histogram viewer"), "#FF00FF"));
      await session.step(136, "And the \"view\" area of line chart viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("line chart viewer"), "#FF00FF"));
      await session.step(137, "And the \"view\" area of bar chart viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("bar chart viewer"), "#FF00FF"));
      await session.step(138, "And the \"view\" area of pie chart viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("pie chart viewer"), "#FF00FF"));
      await session.step(139, "And the \"view\" area of trellis plot viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("trellis plot viewer"), "#FF00FF"));
      await session.step(140, "And the \"view\" area of box plot viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("box plot viewer"), "#FF00FF"));
      await session.step(141, "And the \"view\" area of pie chart viewer should not contain the color \"#2CA02C\"", () => areaNotColor(page, "view", el("pie chart viewer"), "#2CA02C"));
      await session.step(142, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The picker's Cancel puts the color back everywhere", async () => {
      await session.step(145, "When user hovers over \"Asian\" legend item in legend of bar chart viewer", () => hoverOver(page, el("\"Asian\" legend item in legend of bar chart viewer")));
      await session.step(146, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
      await session.step(147, "Then \"Asian\" dialog should be visible", () => shouldBe(page, el("\"Asian\" dialog"), "visible"));
      await session.step(148, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(149, "Then the \"Asian\" item in the legend of pie chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#9467BD"));
      await session.step(150, "When user clicks on CANCEL button in \"Asian\" dialog", () => clickOn(page, el("CANCEL button in \"Asian\" dialog")));
      await session.step(151, "Then \"Asian\" dialog should be absent", () => shouldBe(page, el("\"Asian\" dialog"), "absent"));
      await session.step(152, "And the categorical color of \"Asian\" in \"RACE\" column should be \"#1F77B4\"", () => categoricalColorIs(page, "Asian", "RACE", "#1F77B4"));
      await session.step(153, "And the \"Asian\" item in the legend of bar chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian", el("bar chart viewer"), "#1F77B4"));
      await session.step(154, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#1F77B4"));
      await session.step(155, "And the \"Asian\" item in the legend of scatter plot viewer should be colored \"#1F77B4\"", () => legendItemColor(page, "Asian", el("scatter plot viewer"), "#1F77B4"));
      await session.step(156, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A color picked in the bar chart legend reaches the grid and the other six legends", async () => {
      await session.step(159, "When user hovers over \"Asian\" legend item in legend of bar chart viewer", () => hoverOver(page, el("\"Asian\" legend item in legend of bar chart viewer")));
      await session.step(160, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
      await session.step(161, "Then \"Asian\" dialog should be visible", () => shouldBe(page, el("\"Asian\" dialog"), "visible"));
      await session.step(162, "When user picks the color \"#9467BD\" in the color picker dialog", () => pickColorSwatch(page, "#9467BD"));
      await session.step(163, "And user clicks on OK button in \"Asian\" dialog", () => clickOn(page, el("OK button in \"Asian\" dialog")));
      await session.step(164, "Then \"Asian\" dialog should be absent", () => shouldBe(page, el("\"Asian\" dialog"), "absent"));
      await session.step(165, "And the categorical color of \"Asian\" in \"RACE\" column should be \"#9467BD\"", () => categoricalColorIs(page, "Asian", "RACE", "#9467BD"));
      await session.step(166, "And \"RACE\" column should be color-coded categorically", () => colorCodedCategorically(page, "RACE"));
      await session.step(167, "And the \"color of cell 10 of RACE\" reading of grid should be \"#9467bd\"", () => readingReads(page, "color of cell 10 of RACE", el("grid"), "#9467bd"));
      await session.step(168, "And the \"Asian\" item in the legend of bar chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("bar chart viewer"), "#9467BD"));
      await session.step(169, "And the \"Asian\" item in the legend of scatter plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("scatter plot viewer"), "#9467BD"));
      await session.step(170, "And the \"Asian\" item in the legend of histogram viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("histogram viewer"), "#9467BD"));
      await session.step(171, "And the \"Asian\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("line chart viewer"), "#9467BD"));
      await session.step(172, "And the \"Asian\" item in the legend of pie chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("pie chart viewer"), "#9467BD"));
      await session.step(173, "And the \"Asian\" item in the legend of trellis plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("trellis plot viewer"), "#9467BD"));
      await session.step(174, "And the \"Asian\" item in the legend of box plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("box plot viewer"), "#9467BD"));
      await session.step(175, "And the \"Caucasian\" item in the legend of pie chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("pie chart viewer"), "#FF00FF"));
      await session.step(176, "And the \"pie\" area of pie chart viewer should contain the color \"#9467BD\"", () => areaColor(page, "pie", el("pie chart viewer"), "#9467BD"));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty category gets a \"(no value)\" item whose color can be picked", async () => {
      await session.step(180, "When user adds a calculated column \"RACE under 61\" with formula \"if(${AGE} > 60, null, ${RACE})\"", () => addCalculated(page, "RACE under 61", "if(${AGE} > 60, null, ${RACE})"));
      await session.step(181, "And user sets \"colorColumnName\" property of scatter plot viewer to \"RACE under 61\"", () => setProperty(page, "colorColumnName", el("scatter plot viewer"), "RACE under 61"));
      await session.step(182, "Then the legend of scatter plot viewer should list 5 items", () => legendLists(page, el("scatter plot viewer"), 5));
      await session.step(183, "And \"(no value)\" legend item in legend of scatter plot viewer should be visible", () => shouldBe(page, el("\"(no value)\" legend item in legend of scatter plot viewer"), "visible"));
      await session.step(184, "When user hovers over \"(no value)\" legend item in legend of scatter plot viewer", () => hoverOver(page, el("\"(no value)\" legend item in legend of scatter plot viewer")));
      await session.step(185, "And user clicks on color picker icon", () => clickOn(page, el("color picker icon")));
      await session.step(186, "Then \"(no value)\" dialog should be visible", () => shouldBe(page, el("\"(no value)\" dialog"), "visible"));
      await session.step(187, "And the \"view\" area of scatter plot viewer should not contain the color \"#E377C2\"", () => areaNotColor(page, "view", el("scatter plot viewer"), "#E377C2"));
      await session.step(188, "When user picks the color \"#E377C2\" in the color picker dialog", () => pickColorSwatch(page, "#E377C2"));
      await session.step(189, "And user clicks on OK button in \"(no value)\" dialog", () => clickOn(page, el("OK button in \"(no value)\" dialog")));
      await session.step(190, "Then the \"(no value)\" item in the legend of scatter plot viewer should be colored \"#E377C2\"", () => legendItemColor(page, "(no value)", el("scatter plot viewer"), "#E377C2"));
      await session.step(191, "And the \"view\" area of scatter plot viewer should contain the color \"#E377C2\"", () => areaColor(page, "view", el("scatter plot viewer"), "#E377C2"));
      await session.step(192, "When user sets \"colorColumnName\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
      await session.step(193, "And user removes \"RACE under 61\" column", () => removeColumn(page, "RACE under 61"));
      await session.step(194, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(195, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The palette, the legend column and the visibility come back from a saved layout", async () => {
      await session.step(198, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(199, "And user colors \"RACE\" column categorically:", () => colorCategorical(page, "RACE", [["Caucasian","#2CA02C"],["Asian","#1F77B4"]]), [["Caucasian","#2CA02C"],["Asian","#1F77B4"]]);
      await session.step(202, "And user sets \"colorColumnName\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "colorColumnName", el("scatter plot viewer"), "SEX"));
      await session.step(203, "And user sets \"Legend Visibility\" property of pie chart viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("pie chart viewer"), "Never"));
      await session.step(204, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(205, "And legend of pie chart viewer should be hidden", () => shouldBe(page, el("legend of pie chart viewer"), "hidden"));
      await session.step(206, "And the \"Caucasian\" item in the legend of box plot viewer should be colored \"#2CA02C\"", () => legendItemColor(page, "Caucasian", el("box plot viewer"), "#2CA02C"));
      await session.step(207, "When user loads the saved layout", () => loadLayout(page));
      await session.step(208, "Then \"colorColumnName\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
      await session.step(209, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(210, "And legend of pie chart viewer should be visible", () => shouldBe(page, el("legend of pie chart viewer"), "visible"));
      await session.step(211, "And \"Legend Visibility\" property of pie chart viewer should be \"Always\"", () => propertyShouldBe(page, "Legend Visibility", el("pie chart viewer"), "Always"));
      await session.step(212, "And the categorical color of \"Caucasian\" in \"RACE\" column should be \"#FF00FF\"", () => categoricalColorIs(page, "Caucasian", "RACE", "#FF00FF"));
      await session.step(213, "And the categorical color of \"Asian\" in \"RACE\" column should be \"#9467BD\"", () => categoricalColorIs(page, "Asian", "RACE", "#9467BD"));
      await session.step(214, "And the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF00FF"));
      await session.step(215, "And the \"Caucasian\" item in the legend of histogram viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("histogram viewer"), "#FF00FF"));
      await session.step(216, "And the \"Caucasian\" item in the legend of line chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("line chart viewer"), "#FF00FF"));
      await session.step(217, "And the \"Caucasian\" item in the legend of bar chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("bar chart viewer"), "#FF00FF"));
      await session.step(218, "And the \"Caucasian\" item in the legend of pie chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("pie chart viewer"), "#FF00FF"));
      await session.step(219, "And the \"Caucasian\" item in the legend of trellis plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("trellis plot viewer"), "#FF00FF"));
      await session.step(220, "And the \"Caucasian\" item in the legend of box plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("box plot viewer"), "#FF00FF"));
      await session.step(221, "And the \"Asian\" item in the legend of box plot viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("box plot viewer"), "#9467BD"));
      await session.step(222, "And the \"view\" area of box plot viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("box plot viewer"), "#FF00FF"));
      await session.step(223, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The palette comes back on every viewer from a saved project", async () => {
      await session.step(226, "When user saves the current view as project \"bdd-legend-across-viewers\"", () => saveAsProject(page, "bdd-legend-across-viewers"));
      await session.step(227, "And user closes all views", () => closeAllViews(page));
      await session.step(228, "And user opens the \"bdd-legend-across-viewers\" project", () => openProject(page, "bdd-legend-across-viewers"));
      await session.step(229, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
      await session.step(230, "And the open tableview should have 1 box plot viewer", () => viewerCount(page, 1, "box plot"));
      await session.step(231, "And the categorical color of \"Caucasian\" in \"RACE\" column should be \"#FF00FF\"", () => categoricalColorIs(page, "Caucasian", "RACE", "#FF00FF"));
      await session.step(232, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(233, "And the legend of trellis plot viewer should list 4 items", () => legendLists(page, el("trellis plot viewer"), 4));
      await session.step(234, "And the \"Caucasian\" item in the legend of scatter plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("scatter plot viewer"), "#FF00FF"));
      await session.step(235, "And the \"Caucasian\" item in the legend of histogram viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("histogram viewer"), "#FF00FF"));
      await session.step(236, "And the \"Caucasian\" item in the legend of line chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("line chart viewer"), "#FF00FF"));
      await session.step(237, "And the \"Caucasian\" item in the legend of bar chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("bar chart viewer"), "#FF00FF"));
      await session.step(238, "And the \"Caucasian\" item in the legend of pie chart viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("pie chart viewer"), "#FF00FF"));
      await session.step(239, "And the \"Caucasian\" item in the legend of trellis plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("trellis plot viewer"), "#FF00FF"));
      await session.step(240, "And the \"Caucasian\" item in the legend of box plot viewer should be colored \"#FF00FF\"", () => legendItemColor(page, "Caucasian", el("box plot viewer"), "#FF00FF"));
      await session.step(241, "And the \"Asian\" item in the legend of line chart viewer should be colored \"#9467BD\"", () => legendItemColor(page, "Asian", el("line chart viewer"), "#9467BD"));
      await session.step(242, "And the \"view\" area of pie chart viewer should contain the color \"#FF00FF\"", () => areaColor(page, "view", el("pie chart viewer"), "#FF00FF"));
      await session.step(243, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The coloring goes back where it was found", async () => {
      await session.step(246, "When user removes the coloring of \"RACE\" column", () => colorOff(page, "RACE"));
      await session.step(247, "Then \"RACE\" column should have no color coding", () => noColorCoding(page, "RACE"));
      await session.step(248, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
