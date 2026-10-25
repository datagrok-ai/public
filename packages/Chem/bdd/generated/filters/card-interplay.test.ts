/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/card-interplay.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters]
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
  test("The substructure card with the sketcher, other cards and a reopened panel", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the \"chem-filter\" user settings are put back at feature end", () => userSettingsPutBack(page, "chem-filter"));
    await session.step(19, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("The Structure header dropped on the panel adds an empty card on top", async () => {
      await session.step(24, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(25, "And user adds a card for \"Series\" to the filter panel", () => addCardFor(page, "Series"));
      await session.step(26, "And user drags the \"header Structure\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header Structure", el("grid"), "view", el("filter panel")));
      await session.step(27, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(28, "And the \"cards\" reading of filter panel should be \"Structure, Series\"", () => readingReads(page, "cards", el("filter panel"), "Structure, Series"));
      await session.step(29, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(30, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(31, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared the sketch reaches the card on OK", async () => {
      await session.step(35, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(36, "And user unchecks \"Filter as you draw\" input in sketcher dialog", () => uncheckFilterAsYouDraw(page));
      await session.step(37, "Then \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(38, "When user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(39, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(40, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(41, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(42, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(43, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(44, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(45, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared a further edit waits for OK", async () => {
      await session.step(49, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(50, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(51, "And \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(52, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(53, "And user types \"C1CCNCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNCC1", el("molecule input of sketcher dialog")));
      await session.step(54, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(55, "Then the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(56, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(57, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(58, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(59, "And the filter should pass exactly the molecules of \"Structure\" column containing \"C1CCNCC1\"", () => filterPassesMatching(page, "Structure", "C1CCNCC1"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw checked the edit reaches the grid with the sketcher open", async () => {
      await session.step(63, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(64, "And user checks \"Filter as you draw\" input in sketcher dialog", () => checkFilterAsYouDraw(page));
      await session.step(65, "Then \"Filter as you draw\" input in sketcher dialog should be checked", () => shouldBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(66, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(67, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(68, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(69, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(70, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(71, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(72, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Align and Highlight turned on with the sketcher open repaint the grid and keep the rows", async () => {
      await session.step(76, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(77, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(78, "When user unchecks \"Align\" input in sketcher dialog", () => uncheck(page, el("\"Align\" input in sketcher dialog")));
      await session.step(79, "And user unchecks \"Highlight\" input in sketcher dialog", () => uncheck(page, el("\"Highlight\" input in sketcher dialog")));
      await session.step(80, "Then \"Align\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Align\" input in sketcher dialog"), "checked"));
      await session.step(81, "And \"Highlight\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Highlight\" input in sketcher dialog"), "checked"));
      await session.step(82, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(83, "When user takes a snapshot of grid", () => takeSnapshot(page, el("grid")));
      await session.step(84, "Then grid should not have repainted", () => notRepainted(page, el("grid")));
      await session.step(85, "When user checks \"Align\" input in sketcher dialog", () => check(page, el("\"Align\" input in sketcher dialog")));
      await session.step(86, "And user checks \"Highlight\" input in sketcher dialog", () => check(page, el("\"Highlight\" input in sketcher dialog")));
      await session.step(87, "Then grid should have repainted", () => repainted(page, el("grid")));
      await session.step(88, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(89, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(90, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(91, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(92, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The substructure and category cards intersect and switch off on their own", async () => {
      await session.step(95, "When user clicks on the \"category Pyrrolidines of Series\" area of filter panel", () => clickArea(page, "category Pyrrolidines of Series", el("filter panel")));
      await session.step(96, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(97, "And no rows where \"Series\" is \"Triazoles\" should pass the filter", () => noneOfFiltered(page, "Series", "Triazoles"));
      await session.step(98, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(99, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(100, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(101, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(102, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(103, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(104, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(105, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A switched-off card survives closing and reopening the panel", async () => {
      await session.step(108, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(109, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(110, "Then 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(111, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(112, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(113, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(114, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(115, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(116, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(117, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(118, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(119, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
