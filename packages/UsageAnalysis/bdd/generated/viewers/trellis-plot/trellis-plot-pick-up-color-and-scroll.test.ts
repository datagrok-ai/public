/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-pick-up-color-and-scroll.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, dragAreaBy, hasArea, hasNoArea, menuLists, noErrors, pickFromOpenMenu, propertiesShouldBe, propertyShouldBe, propertyShouldNotBe, readingIs, readingNotAsRemembered, readingReads, rememberReading, reportsNoError, rightClickArea, setProperties, setProperty, viewerCount, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, innerPropertyShouldBe, pickFromViewerMenu, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot Pick Up / Apply, inner color coding and category scrolling", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-pick-up-color-and-scroll.feature", import.meta.url);
  test("Pick Up / Apply carries the Y split, the inner viewer, the legend and the title", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Given user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
    await session.step(20, "And user adds a trellis plot viewer", () => addViewer(page, "trellis plot"));
    await session.step(21, "Then the open tableview should have 2 trellis plot viewers", () => viewerCount(page, 2, "trellis plot"));
    await session.step(22, "When user sets properties of first trellis plot viewer:", () => setProperties(page, el("first trellis plot viewer"), [["Y Column Names","DIS_POP"],["Viewer Type","Bar chart"],["Legend Visibility","Always"],["Legend Position","Left"],["Show Title","true"],["Title","First"]]), [["Y Column Names","DIS_POP"],["Viewer Type","Bar chart"],["Legend Visibility","Always"],["Legend Position","Left"],["Show Title","true"],["Title","First"]]);
    await session.step(29, "Then \"Y Column Names\" property of second trellis plot viewer should not be \"DIS_POP\"", () => propertyShouldNotBe(page, "Y Column Names", el("second trellis plot viewer"), "DIS_POP"));
    await session.step(30, "And second trellis plot viewer should not have a \"y label RA\" area", () => hasNoArea(page, el("second trellis plot viewer"), "y label RA"));
    await session.step(31, "When user picks \"Pick Up / Apply > Pick Up\" from the viewer menu of first trellis plot viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Pick Up", el("first trellis plot viewer")));
    await session.step(32, "And user picks \"Pick Up / Apply > Apply\" from the viewer menu of second trellis plot viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Apply", el("second trellis plot viewer")));
    await session.step(33, "Then properties of second trellis plot viewer should be:", () => propertiesShouldBe(page, el("second trellis plot viewer"), [["Y Column Names","DIS_POP"],["Viewer Type","Bar chart"],["Legend Position","Left"],["Title","First"]]), [["Y Column Names","DIS_POP"],["Viewer Type","Bar chart"],["Legend Position","Left"],["Title","First"]]);
    await session.step(38, "And the \"inner viewer type\" reading of second trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("second trellis plot viewer"), "Bar chart"));
    await session.step(39, "And second trellis plot viewer should have a \"y label RA\" area", () => hasArea(page, el("second trellis plot viewer"), "y label RA"));
    await session.step(40, "And title of second trellis plot viewer should have text \"First\"", () => shouldHaveText(page, el("title of second trellis plot viewer"), "First"));
    await session.step(41, "When user sets \"Y Column Names\" property of first trellis plot viewer to \"RACE\"", () => setProperty(page, "Y Column Names", el("first trellis plot viewer"), "RACE"));
    await session.step(42, "Then the \"y categories\" reading of first trellis plot viewer should be 4", () => readingIs(page, "y categories", el("first trellis plot viewer"), 4));
    await session.step(43, "And the \"y categories\" reading of second trellis plot viewer should be 6", () => readingIs(page, "y categories", el("second trellis plot viewer"), 6));
    await session.step(44, "And \"Y Column Names\" property of second trellis plot viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Y Column Names", el("second trellis plot viewer"), "DIS_POP"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Pie chart inside, colored by RACE, repaints every cell [type=Pie chart, color property=categoryColumnName]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(48, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Pie chart"]]);
    await session.step(52, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
    await session.step(53, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(54, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
    await session.step(55, "And user remembers the \"cell signature F | Black\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Black", el("trellis plot viewer")));
    await session.step(56, "And user remembers the \"cell signature M | Other\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Other", el("trellis plot viewer")));
    await session.step(57, "And user sets \"categoryColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "categoryColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(58, "Then \"categoryColumnName\" inner property of trellis plot viewer should be \"RACE\"", () => innerPropertyShouldBe(page, "categoryColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(59, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(60, "And the \"cell signature M | Asian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
    await session.step(61, "And the \"cell signature F | Black\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Black", el("trellis plot viewer")));
    await session.step(62, "And the \"cell signature M | Other\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Other", el("trellis plot viewer")));
    await session.step(63, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("A Box plot inside, colored by RACE, repaints every cell [type=Box plot, color property=markerColorColumnName]", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(48, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Box plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Box plot"]]);
    await session.step(52, "Then the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
    await session.step(53, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(54, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
    await session.step(55, "And user remembers the \"cell signature F | Black\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Black", el("trellis plot viewer")));
    await session.step(56, "And user remembers the \"cell signature M | Other\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Other", el("trellis plot viewer")));
    await session.step(57, "And user sets \"markerColorColumnName\" inner property of trellis plot viewer to \"RACE\"", () => setInnerProperty(page, "markerColorColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(58, "Then \"markerColorColumnName\" inner property of trellis plot viewer should be \"RACE\"", () => innerPropertyShouldBe(page, "markerColorColumnName", el("trellis plot viewer"), "RACE"));
    await session.step(59, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(60, "And the \"cell signature M | Asian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
    await session.step(61, "And the \"cell signature F | Black\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Black", el("trellis plot viewer")));
    await session.step(62, "And the \"cell signature M | Other\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Other", el("trellis plot viewer")));
    await session.step(63, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
  });
  test("The X selector's menu resets a paged X split", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(72, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX, DIS_POP"],["Y Column Names","RACE"]]), [["X Column Names","SEX, DIS_POP"],["Y Column Names","RACE"]]);
    await session.step(75, "Then the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
    await session.step(76, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
    await session.step(77, "Then the cells of trellis plot viewer should be 6 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 6, 4));
    await session.step(78, "When user right-clicks on the \"x selector 1\" area of trellis plot viewer", () => rightClickArea(page, "x selector 1", el("trellis plot viewer")));
    await session.step(79, "Then the open menu should list \"Reset X columns\"", () => menuLists(page, "Reset X columns"));
    await session.step(80, "When user picks \"Reset X columns\" from the open menu", () => pickFromOpenMenu(page, "Reset X columns"));
    await session.step(81, "Then \"X Column Names\" property of trellis plot viewer should be \"\"", () => propertyShouldBe(page, "X Column Names", el("trellis plot viewer"), ""));
    await session.step(82, "And the \"cells drawn\" reading of trellis plot viewer should be 4", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 4));
    await session.step(83, "And trellis plot viewer should report no error", () => reportsNoError(page, el("trellis plot viewer")));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
  });
  test("Dragging the category scroll handle and the wheel bring other categories in", {tag: ["@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(87, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX, DIS_POP"],["Y Column Names","DIS_POP, RACE"]]), [["X Column Names","SEX, DIS_POP"],["Y Column Names","DIS_POP, RACE"]]);
    await session.step(90, "Then the cells of trellis plot viewer should be 5 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 5));
    await session.step(91, "And trellis plot viewer should have an \"x label F\" area", () => hasArea(page, el("trellis plot viewer"), "x label F"));
    await session.step(92, "And trellis plot viewer should not have an \"x label M\" area", () => hasNoArea(page, el("trellis plot viewer"), "x label M"));
    await session.step(93, "When user drags the \"x scroll handle\" area of trellis plot viewer by 200 pixels to the right", () => dragAreaBy(page, "x scroll handle", el("trellis plot viewer"), 200, "right"));
    await session.step(94, "Then trellis plot viewer should have an \"x label M\" area", () => hasArea(page, el("trellis plot viewer"), "x label M"));
    await session.step(95, "And the cells of trellis plot viewer should be 5 wide and 5 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 5));
    await session.step(96, "And trellis plot viewer should not have a \"y label PsA\" area", () => hasNoArea(page, el("trellis plot viewer"), "y label PsA"));
    await session.step(97, "When user scrolls the mouse wheel down 3 times over the \"view\" area of trellis plot viewer", () => wheelOverAreaTimes(page, "down", 3, "view", el("trellis plot viewer")));
    await session.step(98, "Then trellis plot viewer should have a \"y label PsA\" area", () => hasArea(page, el("trellis plot viewer"), "y label PsA"));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
  });
});
