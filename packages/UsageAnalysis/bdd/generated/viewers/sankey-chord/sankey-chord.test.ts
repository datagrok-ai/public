/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sankey-chord/sankey-chord.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sankey, charts.viewer.chord]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {boundTable, clickArea, hasArea, hasNoArea, hoverArea, noErrors, pointerAway, propertiesShouldBe, propertyShouldBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {offersColumns, readingContains, readingNotContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sankey and Chord columns, and redrawing on every filter change", () => {
  const session = feature(test, "features/viewers/sankey-chord/sankey-chord.feature", import.meta.url);
  test("Sankey starts on SEX, RACE and AGE and takes other columns (GROK-18048)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(20, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(21, "When user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(22, "Then sankey viewer should be visible", () => shouldBe(page, el("sankey viewer"), "visible"));
    await session.step(23, "And sankey viewer should be bound to table \"demog\"", () => boundTable(page, el("sankey viewer"), "demog"));
    await session.step(24, "And the \"node names\" reading of sankey viewer should be \"F, M, Caucasian, Other, Asian, Black\"", () => readingReads(page, "node names", el("sankey viewer"), "F, M, Caucasian, Other, Asian, Black"));
    await session.step(25, "And the \"links\" reading of sankey viewer should be 5850", () => readingIs(page, "links", el("sankey viewer"), 5850));
    await session.step(26, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(27, "And user clicks on settings icon of sankey viewer", () => clickOn(page, el("settings icon of sankey viewer")));
    await session.step(28, "Then \"Source\" property in context panel should be visible", () => shouldBe(page, el("\"Source\" property in context panel"), "visible"));
    await session.step(29, "And properties of sankey viewer should be:", () => propertiesShouldBe(page, el("sankey viewer"), [["Source","SEX"],["Target","RACE"],["Value","AGE"]]), [["Source","SEX"],["Target","RACE"],["Value","AGE"]]);
    await session.step(33, "And \"Source\" property in context panel should offer the columns \"USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY\"", () => offersColumns(page, el("\"Source\" property in context panel"), "USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY"));
    await session.step(34, "And \"Target\" property in context panel should offer the columns \"USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY\"", () => offersColumns(page, el("\"Target\" property in context panel"), "USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY"));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await session.step(36, "When user selects \"RACE\" in \"Source\" property in context panel", () => selectIn(page, "RACE", el("\"Source\" property in context panel")));
    await session.step(37, "Then \"Source\" property of sankey viewer should be \"RACE\"", () => propertyShouldBe(page, "Source", el("sankey viewer"), "RACE"));
    await session.step(38, "When user selects \"DIS_POP\" in \"Target\" property in context panel", () => selectIn(page, "DIS_POP", el("\"Target\" property in context panel")));
    await session.step(39, "Then \"Target\" property of sankey viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Target", el("sankey viewer"), "DIS_POP"));
    await session.step(40, "When user selects \"WEIGHT\" in \"Value\" property in context panel", () => selectIn(page, "WEIGHT", el("\"Value\" property in context panel")));
    await session.step(41, "Then properties of sankey viewer should be:", () => propertiesShouldBe(page, el("sankey viewer"), [["Source","RACE"],["Target","DIS_POP"],["Value","WEIGHT"]]), [["Source","RACE"],["Target","DIS_POP"],["Value","WEIGHT"]]);
    await session.step(45, "And the \"node names\" reading of sankey viewer should contain \"Psoriasis\"", () => readingContains(page, "node names", el("sankey viewer"), "Psoriasis"));
    await session.step(46, "And the \"node names\" reading of sankey viewer should contain \"Caucasian\"", () => readingContains(page, "node names", el("sankey viewer"), "Caucasian"));
    await session.step(47, "And the \"node names\" reading of sankey viewer should not contain \"M\"", () => readingNotContains(page, "node names", el("sankey viewer"), "M"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sankey follows the filter (GROK-18035)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(51, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(52, "And user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(53, "Then the \"node names\" reading of sankey viewer should contain \"M\"", () => readingContains(page, "node names", el("sankey viewer"), "M"));
    await session.step(54, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(55, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(56, "When user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
    await session.step(57, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
    await session.step(58, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
    await session.step(59, "And the \"links\" reading of sankey viewer should be 3243", () => readingIs(page, "links", el("sankey viewer"), 3243));
    await session.step(60, "And the \"node names\" reading of sankey viewer should contain \"F\"", () => readingContains(page, "node names", el("sankey viewer"), "F"));
    await session.step(61, "And the \"node names\" reading of sankey viewer should not contain \"M\"", () => readingNotContains(page, "node names", el("sankey viewer"), "M"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
    await session.step(63, "When user hovers over the \"link F -> Caucasian\" area of sankey viewer", () => hoverArea(page, "link F -> Caucasian", el("sankey viewer")));
    await session.step(64, "Then tooltip should contain text \"2823 rows\"", () => shouldContainText(page, el("tooltip"), "2823 rows"));
    await session.step(65, "When user hovers over the \"link F -> Asian\" area of sankey viewer", () => hoverArea(page, "link F -> Asian", el("sankey viewer")));
    await session.step(66, "Then tooltip should contain text \"37 rows\"", () => shouldContainText(page, el("tooltip"), "37 rows"));
    await session.step(67, "When user moves the pointer away from sankey viewer", () => pointerAway(page, el("sankey viewer")));
    await session.step(68, "Then no errors should have been logged", () => noErrors(page));
    await session.step(69, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(70, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(71, "Then 5850 rows should pass the filter", () => filterPasses(page, 5850));
    await session.step(72, "And the \"links\" reading of sankey viewer should be 5850", () => readingIs(page, "links", el("sankey viewer"), 5850));
    await session.step(73, "And the \"node names\" reading of sankey viewer should contain \"M\"", () => readingContains(page, "node names", el("sankey viewer"), "M"));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sankey with a filter no row passes draws nothing and logs no error (GROK-21110)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(77, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(78, "And user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(79, "Then the \"links\" reading of sankey viewer should be 5850", () => readingIs(page, "links", el("sankey viewer"), 5850));
    await session.step(80, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(81, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(82, "When user clicks on the \"category true of CONTROL\" area of filter panel", () => clickArea(page, "category true of CONTROL", el("filter panel")));
    await session.step(83, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(84, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(85, "And the \"links\" reading of sankey viewer should be 0", () => readingIs(page, "links", el("sankey viewer"), 0));
    await session.step(86, "And the \"nodes\" reading of sankey viewer should be 0", () => readingIs(page, "nodes", el("sankey viewer"), 0));
    await session.step(87, "And no errors should have been logged", () => noErrors(page));
  });
  test("Chord takes From set to the column To holds, and redraws on a filter change without a click (GROK-21111, GROK-17772)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(90, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(91, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(92, "When user clicks on first \"Chord\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Chord\" card in \"Add Viewer\" dialog")));
    await session.step(93, "Then chord viewer should be visible", () => shouldBe(page, el("chord viewer"), "visible"));
    await session.step(94, "And the \"categories\" reading of chord viewer should be 6", () => readingIs(page, "categories", el("chord viewer"), 6));
    await session.step(95, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(96, "And user clicks on settings icon of chord viewer", () => clickOn(page, el("settings icon of chord viewer")));
    await session.step(97, "Then \"From\" property in context panel should be visible", () => shouldBe(page, el("\"From\" property in context panel"), "visible"));
    await session.step(98, "And properties of chord viewer should be:", () => propertiesShouldBe(page, el("chord viewer"), [["From","SEX"],["To","RACE"]]), [["From","SEX"],["To","RACE"]]);
    await session.step(101, "And no errors should have been logged", () => noErrors(page));
    await session.step(102, "When user selects \"RACE\" in \"From\" property in context panel", () => selectIn(page, "RACE", el("\"From\" property in context panel")));
    await session.step(103, "Then \"From\" property of chord viewer should be \"RACE\"", () => propertyShouldBe(page, "From", el("chord viewer"), "RACE"));
    await session.step(104, "And the \"categories\" reading of chord viewer should be 4", () => readingIs(page, "categories", el("chord viewer"), 4));
    await session.step(105, "And no errors should have been logged", () => noErrors(page));
    await session.step(106, "When user selects \"DIS_POP\" in \"To\" property in context panel", () => selectIn(page, "DIS_POP", el("\"To\" property in context panel")));
    await session.step(107, "Then properties of chord viewer should be:", () => propertiesShouldBe(page, el("chord viewer"), [["From","RACE"],["To","DIS_POP"]]), [["From","RACE"],["To","DIS_POP"]]);
    await session.step(110, "And chord viewer should have a \"category Black\" area", () => hasArea(page, el("chord viewer"), "category Black"));
    await session.step(111, "And chord viewer should have a \"category PsA\" area", () => hasArea(page, el("chord viewer"), "category PsA"));
    await session.step(112, "And chord viewer should not have a \"category F\" area", () => hasNoArea(page, el("chord viewer"), "category F"));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
    await session.step(114, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(115, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(116, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(117, "Then the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(118, "And chord viewer should have a \"category Asian\" area", () => hasArea(page, el("chord viewer"), "category Asian"));
    await session.step(119, "And chord viewer should not have a \"category Black\" area", () => hasNoArea(page, el("chord viewer"), "category Black"));
    await session.step(120, "And no errors should have been logged", () => noErrors(page));
    await session.step(121, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(122, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(123, "Then 5850 rows should pass the filter", () => filterPasses(page, 5850));
    await session.step(124, "And chord viewer should have a \"category Black\" area", () => hasArea(page, el("chord viewer"), "category Black"));
    await session.step(125, "And no errors should have been logged", () => noErrors(page));
  });
});
