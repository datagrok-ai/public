/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/card-interplay.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters, public#697, public#3328]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {checkFilterAsYouDraw, uncheckFilterAsYouDraw} from '../../bindings/dialogs.js';
import {filterPassesMatching} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, hoverOver, pressKeyIn, shouldBe, shouldNotBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneOfFiltered, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, noErrors, notRepainted, readingReads, repainted, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The substructure card with the sketcher, other cards and a reopened panel", () => {
  const session = feature(test, "features/filters/card-interplay.feature", import.meta.url);
  test("The substructure card with the sketcher, other cards and a reopened panel", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters", "@serial", "@realizes:public#697", "@realizes:public#3328"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the \"chem-filter\" user settings are put back at feature end", () => userSettingsPutBack(page, "chem-filter"));
    await session.step(22, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(23, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(24, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("The Structure header dropped on the panel adds an empty card on top", async () => {
      await session.step(27, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(28, "And user adds a card for \"Series\" to the filter panel", () => addCardFor(page, "Series"));
      await session.step(29, "And user drags the \"header Structure\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header Structure", el("grid"), "view", el("filter panel")));
      await session.step(30, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(31, "And the \"cards\" reading of filter panel should be \"Structure, Series\"", () => readingReads(page, "cards", el("filter panel"), "Structure, Series"));
      await session.step(32, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(33, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(34, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cancel on a new structure drawn with Filter as you draw gives every row back (public#697)", async () => {
      await session.step(39, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(40, "And user checks \"Filter as you draw\" input in sketcher dialog", () => checkFilterAsYouDraw(page));
      await session.step(41, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(42, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(43, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(44, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
      await session.step(45, "Then sketcher dialog should be hidden", () => shouldBe(page, el("sketcher dialog"), "hidden"));
      await session.step(46, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(47, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(48, "And the \"filtering of Structure\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "false"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared the sketch reaches the card on OK", async () => {
      await session.step(52, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(53, "And user unchecks \"Filter as you draw\" input in sketcher dialog", () => uncheckFilterAsYouDraw(page));
      await session.step(54, "Then \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(55, "When user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(56, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(57, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(58, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(59, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(60, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(61, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(62, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared a further edit waits for OK", async () => {
      await session.step(66, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(67, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(68, "And \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(69, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(70, "And user types \"C1CCNCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNCC1", el("molecule input of sketcher dialog")));
      await session.step(71, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(72, "Then the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(73, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(74, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(75, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(76, "And the filter should pass exactly the molecules of \"Structure\" column containing \"C1CCNCC1\"", () => filterPassesMatching(page, "Structure", "C1CCNCC1"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw checked the edit reaches the grid with the sketcher open", async () => {
      await session.step(80, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(81, "And user checks \"Filter as you draw\" input in sketcher dialog", () => checkFilterAsYouDraw(page));
      await session.step(82, "Then \"Filter as you draw\" input in sketcher dialog should be checked", () => shouldBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(83, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(84, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(85, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(86, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(87, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(88, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(89, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(90, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cancel on an edited structure drawn with Filter as you draw gives the old structure back (public#697)", async () => {
      await session.step(94, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(95, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(96, "And \"Filter as you draw\" input in sketcher dialog should be checked", () => shouldBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(97, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(98, "And user types \"C1CCNCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNCC1", el("molecule input of sketcher dialog")));
      await session.step(99, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(100, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(101, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
      await session.step(102, "Then sketcher dialog should be hidden", () => shouldBe(page, el("sketcher dialog"), "hidden"));
      await session.step(103, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(104, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Align and Highlight turned on with the sketcher open repaint the grid and keep the rows", async () => {
      await session.step(108, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(109, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(110, "When user unchecks \"Align\" input in sketcher dialog", () => uncheck(page, el("\"Align\" input in sketcher dialog")));
      await session.step(111, "And user unchecks \"Highlight\" input in sketcher dialog", () => uncheck(page, el("\"Highlight\" input in sketcher dialog")));
      await session.step(112, "Then \"Align\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Align\" input in sketcher dialog"), "checked"));
      await session.step(113, "And \"Highlight\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Highlight\" input in sketcher dialog"), "checked"));
      await session.step(114, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(115, "When user takes a snapshot of grid", () => takeSnapshot(page, el("grid")));
      await session.step(116, "Then grid should not have repainted", () => notRepainted(page, el("grid")));
      await session.step(117, "When user checks \"Align\" input in sketcher dialog", () => check(page, el("\"Align\" input in sketcher dialog")));
      await session.step(118, "And user checks \"Highlight\" input in sketcher dialog", () => check(page, el("\"Highlight\" input in sketcher dialog")));
      await session.step(119, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(120, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(121, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(122, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(123, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(124, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The substructure and category cards intersect and switch off on their own", async () => {
      await session.step(127, "When user clicks on the \"category Pyrrolidines of Series\" area of filter panel", () => clickArea(page, "category Pyrrolidines of Series", el("filter panel")));
      await session.step(128, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(129, "And no rows where \"Series\" is \"Triazoles\" should pass the filter", () => noneOfFiltered(page, "Series", "Triazoles"));
      await session.step(130, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(131, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(132, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(133, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(134, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(135, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(136, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A switched-off card survives closing and reopening the panel", async () => {
      await session.step(140, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(141, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(142, "Then 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(143, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(144, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(145, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(146, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(147, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(148, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(149, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(150, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(151, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(152, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The panel switched off and on again filters by the structure the card holds (public#3328)", async () => {
      await session.step(156, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(157, "And user unchecks master of filter panel", () => uncheck(page, el("master of filter panel")));
      await session.step(158, "Then the \"active\" reading of filter panel should be \"false\"", () => readingReads(page, "active", el("filter panel"), "false"));
      await session.step(159, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(160, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(161, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(162, "And user checks master of filter panel", () => check(page, el("master of filter panel")));
      await session.step(163, "Then the \"active\" reading of filter panel should be \"true\"", () => readingReads(page, "active", el("filter panel"), "true"));
      await session.step(164, "And 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(165, "And the \"filtering of Structure\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of Structure", el("filter panel"), "true"));
      await session.step(166, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
