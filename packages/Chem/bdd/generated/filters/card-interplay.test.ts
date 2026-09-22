/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/card-interplay.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {filterPassesMatching} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, hoverOver, pressKeyIn, shouldBe, shouldNotBe, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneOfFiltered, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The substructure card with the sketcher, other cards and a reopened panel", () => {
  const session = feature(test, "features/filters/card-interplay.feature", import.meta.url);
  test("The substructure card with the sketcher, other cards and a reopened panel", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("The Structure header dropped on the panel adds an empty card on top", async () => {
      await session.step(18, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(19, "And user adds a card for \"Series\" to the filter panel", () => addCardFor(page, "Series"));
      await session.step(20, "And user drags the \"header Structure\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header Structure", el("grid"), "view", el("filter panel")));
      await session.step(21, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(22, "And the \"cards\" reading of filter panel should be \"Structure, Series\"", () => readingReads(page, "cards", el("filter panel"), "Structure, Series"));
      await session.step(23, "And the \"type of Structure\" reading of filter panel should be \"Chem:substructureFilter\"", () => readingReads(page, "type of Structure", el("filter panel"), "Chem:substructureFilter"));
      await session.step(24, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(25, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared the sketch reaches the card on OK", async () => {
      await session.step(29, "When user clicks on \"Sketch\" text in \"Structure\" filter card", () => clickOn(page, el("\"Sketch\" text in \"Structure\" filter card")));
      await session.step(30, "And user unchecks \"Filter as you draw\" input in sketcher dialog", () => uncheck(page, el("\"Filter as you draw\" input in sketcher dialog")));
      await session.step(31, "Then \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(32, "When user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(33, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(34, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(35, "And the \"structure of Structure\" reading of filter panel should be \"\"", () => readingReads(page, "structure of Structure", el("filter panel"), ""));
      await session.step(36, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(37, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(38, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(39, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw cleared a further edit waits for OK", async () => {
      await session.step(43, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(44, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(45, "And \"Filter as you draw\" input in sketcher dialog should not be checked", () => shouldNotBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(46, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(47, "And user types \"C1CCNCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCNCC1", el("molecule input of sketcher dialog")));
      await session.step(48, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(49, "Then the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(50, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(51, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(52, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(53, "And the filter should pass exactly the molecules of \"Structure\" column containing \"C1CCNCC1\"", () => filterPassesMatching(page, "Structure", "C1CCNCC1"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw checked the edit reaches the grid with the sketcher open", async () => {
      await session.step(57, "When user clicks on the \"card Structure\" area of filter panel", () => clickArea(page, "card Structure", el("filter panel")));
      await session.step(58, "And user checks \"Filter as you draw\" input in sketcher dialog", () => check(page, el("\"Filter as you draw\" input in sketcher dialog")));
      await session.step(59, "Then \"Filter as you draw\" input in sketcher dialog should be checked", () => shouldBe(page, el("\"Filter as you draw\" input in sketcher dialog"), "checked"));
      await session.step(60, "When user clears molecule input of sketcher dialog", () => clearField(page, el("molecule input of sketcher dialog")));
      await session.step(61, "And user types \"c1ccncc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccncc1", el("molecule input of sketcher dialog")));
      await session.step(62, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(63, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(64, "And sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(65, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(66, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The substructure and category cards intersect and switch off on their own", async () => {
      await session.step(70, "When user clicks on the \"category Pyrrolidines of Series\" area of filter panel", () => clickArea(page, "category Pyrrolidines of Series", el("filter panel")));
      await session.step(71, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(72, "And no rows where \"Series\" is \"Triazoles\" should pass the filter", () => noneOfFiltered(page, "Series", "Triazoles"));
      await session.step(73, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(74, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(75, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(76, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(77, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(78, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(79, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A switched-off card survives closing and reopening the panel", async () => {
      await session.step(83, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(84, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(85, "Then 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(86, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(87, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(88, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(89, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(90, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(91, "And 21 rows should pass the filter", () => filterPasses(page, 21));
      await session.step(92, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(93, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(94, "Then 5 rows should pass the filter", () => filterPasses(page, 5));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
