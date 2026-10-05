/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sankey-chord/sankey-chord.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sankey, charts.viewer.chord]
--- */
import {test} from '@playwright/test';
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
import {boundTable, clickArea, noErrors, pointerAway, propertiesShouldBe, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Sankey and Chord columns, and redrawing on every filter change", () => {
  const session = feature(test, "features/viewers/sankey-chord/sankey-chord.feature", import.meta.url);
  test("Sankey starts on SEX, RACE and AGE and takes other columns (GROK-18048)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(26, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(27, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(28, "When user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(29, "Then sankey viewer should be visible", () => shouldBe(page, el("sankey viewer"), "visible"));
    await session.step(30, "And sankey viewer should be bound to table \"demog\"", () => boundTable(page, el("sankey viewer"), "demog"));
    await session.step(31, "And \"Caucasian\" text in sankey viewer should be visible", () => shouldBe(page, el("\"Caucasian\" text in sankey viewer"), "visible"));
    await session.step(32, "And \"M\" text in sankey viewer should be visible", () => shouldBe(page, el("\"M\" text in sankey viewer"), "visible"));
    await session.step(33, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(34, "And user clicks on settings icon of sankey viewer", () => clickOn(page, el("settings icon of sankey viewer")));
    await session.step(35, "Then \"Source\" property in context panel should be visible", () => shouldBe(page, el("\"Source\" property in context panel"), "visible"));
    await session.step(36, "And properties of sankey viewer should be:", () => propertiesShouldBe(page, el("sankey viewer"), [["Source","SEX"],["Target","RACE"],["Value","AGE"]]), [["Source","SEX"],["Target","RACE"],["Value","AGE"]]);
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
    await session.step(41, "When user selects \"RACE\" in \"Source\" property in context panel", () => selectIn(page, "RACE", el("\"Source\" property in context panel")));
    await session.step(42, "Then \"Source\" property of sankey viewer should be \"RACE\"", () => propertyShouldBe(page, "Source", el("sankey viewer"), "RACE"));
    await session.step(43, "When user selects \"DIS_POP\" in \"Target\" property in context panel", () => selectIn(page, "DIS_POP", el("\"Target\" property in context panel")));
    await session.step(44, "Then \"Target\" property of sankey viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "Target", el("sankey viewer"), "DIS_POP"));
    await session.step(45, "When user selects \"WEIGHT\" in \"Value\" property in context panel", () => selectIn(page, "WEIGHT", el("\"Value\" property in context panel")));
    await session.step(46, "Then properties of sankey viewer should be:", () => propertiesShouldBe(page, el("sankey viewer"), [["Source","RACE"],["Target","DIS_POP"],["Value","WEIGHT"]]), [["Source","RACE"],["Target","DIS_POP"],["Value","WEIGHT"]]);
    await session.step(50, "And \"Psoriasis\" text in sankey viewer should be visible", () => shouldBe(page, el("\"Psoriasis\" text in sankey viewer"), "visible"));
    await session.step(51, "And \"M\" text in sankey viewer should be absent", () => shouldBe(page, el("\"M\" text in sankey viewer"), "absent"));
    await session.step(52, "And \"Caucasian\" text in sankey viewer should be visible", () => shouldBe(page, el("\"Caucasian\" text in sankey viewer"), "visible"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sankey follows the filter (GROK-18035)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(56, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(57, "And user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(58, "Then \"M\" text in sankey viewer should be visible", () => shouldBe(page, el("\"M\" text in sankey viewer"), "visible"));
    await session.step(59, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(60, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(61, "When user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
    await session.step(62, "Then 3243 rows should pass the filter", () => filterPasses(page, 3243));
    await session.step(63, "And the filter should pass exactly the rows where \"SEX\" is \"F\"", () => filterIsExactlyCategory(page, "SEX", "F"));
    await session.step(64, "And \"M\" text in sankey viewer should be absent", () => shouldBe(page, el("\"M\" text in sankey viewer"), "absent"));
    await session.step(65, "And \"F\" text in sankey viewer should be visible", () => shouldBe(page, el("\"F\" text in sankey viewer"), "visible"));
    await session.step(66, "And \"M\" text in sankey viewer should be absent", () => shouldBe(page, el("\"M\" text in sankey viewer"), "absent"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
    await session.step(68, "When user hovers over sankey viewer", () => hoverOver(page, el("sankey viewer")));
    await session.step(69, "Then tooltip should contain text \"2823 rows\"", () => shouldContainText(page, el("tooltip"), "2823 rows"));
    await session.step(70, "When user moves the pointer away from sankey viewer", () => pointerAway(page, el("sankey viewer")));
    await session.step(71, "Then no errors should have been logged", () => noErrors(page));
    await session.step(72, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(73, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(74, "Then 5850 rows should pass the filter", () => filterPasses(page, 5850));
    await session.step(75, "And \"F\" text in sankey viewer should be visible", () => shouldBe(page, el("\"F\" text in sankey viewer"), "visible"));
    await session.step(76, "And \"M\" text in sankey viewer should be visible", () => shouldBe(page, el("\"M\" text in sankey viewer"), "visible"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
  });
  test("Sankey with a filter no row passes logs no error (GROK-21110)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await knownFailure(async () => {
      await session.step(81, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
      await session.step(82, "And user clicks on first \"Sankey\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
      await session.step(83, "Then \"M\" text in sankey viewer should be visible", () => shouldBe(page, el("\"M\" text in sankey viewer"), "visible"));
      await session.step(84, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(85, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(86, "When user clicks on the \"category true of CONTROL\" area of filter panel", () => clickArea(page, "category true of CONTROL", el("filter panel")));
      await session.step(87, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(88, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(89, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("Chord redraws on a filter change without a click (GROK-17772)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(92, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(93, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(94, "When user clicks on first \"Chord\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Chord\" card in \"Add Viewer\" dialog")));
    await session.step(95, "Then chord viewer should be visible", () => shouldBe(page, el("chord viewer"), "visible"));
    await session.step(96, "And \"Asian\" text in chord viewer should be visible", () => shouldBe(page, el("\"Asian\" text in chord viewer"), "visible"));
    await session.step(97, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(98, "And user clicks on settings icon of chord viewer", () => clickOn(page, el("settings icon of chord viewer")));
    await session.step(99, "Then \"From\" property in context panel should be visible", () => shouldBe(page, el("\"From\" property in context panel"), "visible"));
    await session.step(100, "And properties of chord viewer should be:", () => propertiesShouldBe(page, el("chord viewer"), [["From","SEX"],["To","RACE"]]), [["From","SEX"],["To","RACE"]]);
    await session.step(103, "And no errors should have been logged", () => noErrors(page));
    await session.step(104, "When user selects \"DIS_POP\" in \"To\" property in context panel", () => selectIn(page, "DIS_POP", el("\"To\" property in context panel")));
    await session.step(105, "Then \"To\" property of chord viewer should be \"DIS_POP\"", () => propertyShouldBe(page, "To", el("chord viewer"), "DIS_POP"));
    await session.step(106, "And \"PsA\" text in chord viewer should be visible", () => shouldBe(page, el("\"PsA\" text in chord viewer"), "visible"));
    await session.step(107, "When user selects \"RACE\" in \"From\" property in context panel", () => selectIn(page, "RACE", el("\"From\" property in context panel")));
    await session.step(108, "Then properties of chord viewer should be:", () => propertiesShouldBe(page, el("chord viewer"), [["From","RACE"],["To","DIS_POP"]]), [["From","RACE"],["To","DIS_POP"]]);
    await session.step(111, "And \"F\" text in chord viewer should be absent", () => shouldBe(page, el("\"F\" text in chord viewer"), "absent"));
    await session.step(112, "And \"Black\" text in chord viewer should be visible", () => shouldBe(page, el("\"Black\" text in chord viewer"), "visible"));
    await session.step(113, "And \"PsA\" text in chord viewer should be visible", () => shouldBe(page, el("\"PsA\" text in chord viewer"), "visible"));
    await session.step(114, "And \"F\" text in chord viewer should be absent", () => shouldBe(page, el("\"F\" text in chord viewer"), "absent"));
    await session.step(115, "And no errors should have been logged", () => noErrors(page));
    await session.step(116, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(117, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(118, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(119, "Then the filter should pass exactly the rows where \"RACE\" is \"Asian\"", () => filterIsExactlyCategory(page, "RACE", "Asian"));
    await session.step(120, "And \"Black\" text in chord viewer should be absent", () => shouldBe(page, el("\"Black\" text in chord viewer"), "absent"));
    await session.step(121, "And \"Asian\" text in chord viewer should be visible", () => shouldBe(page, el("\"Asian\" text in chord viewer"), "visible"));
    await session.step(122, "And \"Black\" text in chord viewer should be absent", () => shouldBe(page, el("\"Black\" text in chord viewer"), "absent"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
    await session.step(124, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(125, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
    await session.step(126, "Then 5850 rows should pass the filter", () => filterPasses(page, 5850));
    await session.step(127, "And \"Black\" text in chord viewer should be visible", () => shouldBe(page, el("\"Black\" text in chord viewer"), "visible"));
    await session.step(128, "And no errors should have been logged", () => noErrors(page));
  });
  test("Setting the Chord's From to the column To holds logs no error (GROK-21111)", {tag: ["@viewers", "@realizes:charts.viewer.sankey", "@realizes:charts.viewer.chord", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await knownFailure(async () => {
      await session.step(132, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
      await session.step(133, "And user clicks on first \"Chord\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Chord\" card in \"Add Viewer\" dialog")));
      await session.step(134, "Then chord viewer should be visible", () => shouldBe(page, el("chord viewer"), "visible"));
      await session.step(135, "When user clicks on grid", () => clickOn(page, el("grid")));
      await session.step(136, "And user clicks on settings icon of chord viewer", () => clickOn(page, el("settings icon of chord viewer")));
      await session.step(137, "Then \"From\" property in context panel should be visible", () => shouldBe(page, el("\"From\" property in context panel"), "visible"));
      await session.step(138, "When user selects \"RACE\" in \"From\" property in context panel", () => selectIn(page, "RACE", el("\"From\" property in context panel")));
      await session.step(139, "Then \"From\" property of chord viewer should be \"RACE\"", () => propertyShouldBe(page, "From", el("chord viewer"), "RACE"));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
  });
});
