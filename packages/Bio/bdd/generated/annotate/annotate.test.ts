/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/annotate/annotate.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.analyze.compare-sequences, bio.annotate.scan-liabilities, bio.annotate.manage-annotations]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {carriesAtLeast, carriesNone, countBelowHits, hitsMatch, listMatchesAnnotations, oneFewer} from '../../bindings/annotations.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType, columnTag, columnType, hasNoColumn, joinedValues} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnNamed, newColumnsCount, noNewColumn, pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openDatasetRowsAs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Comparing, scanning and annotating antibody sequences", () => {
  const session = feature(test, "features/annotate/annotate.feature", import.meta.url);
  test("Comparing, scanning and annotating antibody sequences", {tag: ["@journey", "@realizes:bio.analyze.compare-sequences", "@realizes:bio.annotate.scan-liabilities", "@realizes:bio.annotate.manage-annotations"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(8, "Given user is logged in", () => loggedIn(page));
    await session.step(9, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(10, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(11, "Then \"AntibodyHC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyHC", "Macromolecule"));
    await session.step(12, "And \"AntibodyLC\" column should have semantic type \"Macromolecule\"", () => columnSemType(page, "AntibodyLC", "Macromolecule"));
    await run.scenario("Compare sequences pairs the two chains with the defaults", async () => {
      await session.step(15, "When user picks \"Bio > Analyze > Compare sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Compare sequences..."));
      await session.step(16, "Then \"Compare Sequences\" dialog should be visible", () => shouldBe(page, el("\"Compare Sequences\" dialog"), "visible"));
      await session.step(17, "And \"Sequence column 1\" input in \"Compare Sequences\" dialog should have value \"AntibodyHC\"", () => shouldHaveValue(page, el("\"Sequence column 1\" input in \"Compare Sequences\" dialog"), "AntibodyHC"));
      await session.step(18, "And \"Sequence column 2\" input in \"Compare Sequences\" dialog should have value \"AntibodyLC\"", () => shouldHaveValue(page, el("\"Sequence column 2\" input in \"Compare Sequences\" dialog"), "AntibodyLC"));
      await session.step(19, "When user clicks on OK button in \"Compare Sequences\" dialog", () => clickOn(page, el("OK button in \"Compare Sequences\" dialog")));
      await session.step(20, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(21, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(22, "And a new column \"AntibodyHC vs AntibodyLC\" should have been added", () => newColumnNamed(page, "AntibodyHC vs AntibodyLC"));
      await session.step(23, "And \"AntibodyHC vs AntibodyLC\" column should have tag \"cell.renderer\" equal to \"MacromoleculeDifference\"", () => columnTag(page, "AntibodyHC vs AntibodyLC", "cell.renderer", "MacromoleculeDifference"));
      await session.step(24, "And every value of \"AntibodyHC vs AntibodyLC\" column should be \"AntibodyHC\" and \"AntibodyLC\" of the same row joined by \"#\"", () => joinedValues(page, "AntibodyHC vs AntibodyLC", "AntibodyHC", "AntibodyLC", "#"));
      await session.step(25, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A custom name and swapped columns are honoured", async () => {
      await session.step(28, "When user picks \"Bio > Analyze > Compare sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Compare sequences..."));
      await session.step(29, "And user selects \"AntibodyLC\" in \"Sequence column 1\" input in \"Compare Sequences\" dialog", () => selectIn(page, "AntibodyLC", el("\"Sequence column 1\" input in \"Compare Sequences\" dialog")));
      await session.step(30, "And user selects \"AntibodyHC\" in \"Sequence column 2\" input in \"Compare Sequences\" dialog", () => selectIn(page, "AntibodyHC", el("\"Sequence column 2\" input in \"Compare Sequences\" dialog")));
      await session.step(31, "And user enters \"LC minus HC\" into \"Result column name\" input in \"Compare Sequences\" dialog", () => enterInto(page, "LC minus HC", el("\"Result column name\" input in \"Compare Sequences\" dialog")));
      await session.step(32, "And user clicks on OK button in \"Compare Sequences\" dialog", () => clickOn(page, el("OK button in \"Compare Sequences\" dialog")));
      await session.step(33, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(34, "And 1 new column should have been added", () => newColumnsCount(page, 1));
      await session.step(35, "And a new column \"LC minus HC\" should have been added", () => newColumnNamed(page, "LC minus HC"));
      await session.step(36, "And every value of \"LC minus HC\" column should be \"AntibodyLC\" and \"AntibodyHC\" of the same row joined by \"#\"", () => joinedValues(page, "LC minus HC", "AntibodyLC", "AntibodyHC", "#"));
    });
    await run.scenario("The same column twice is refused and adds nothing", async () => {
      await session.step(39, "When user picks \"Bio > Analyze > Compare sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Compare sequences..."));
      await session.step(40, "And user selects \"AntibodyHC\" in \"Sequence column 2\" input in \"Compare Sequences\" dialog", () => selectIn(page, "AntibodyHC", el("\"Sequence column 2\" input in \"Compare Sequences\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Compare Sequences\" dialog", () => clickOn(page, el("OK button in \"Compare Sequences\" dialog")));
      await session.step(42, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(43, "And an error balloon containing \"distinct columns\" should have been shown", () => errorBalloonText(page, "distinct columns"));
      await session.step(44, "And no new column should have been added", () => noNewColumn(page));
    });
    await run.scenario("Scan Liabilities opens with its documented rule defaults", async () => {
      await session.step(47, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
      await session.step(48, "Then \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "visible"));
      await session.step(49, "And \"Deamidation (NG)\" checkbox in \"Scan Sequence Liabilities\" dialog should be checked", () => shouldBe(page, el("\"Deamidation (NG)\" checkbox in \"Scan Sequence Liabilities\" dialog"), "checked"));
      await session.step(50, "And \"Free Cysteine\" checkbox in \"Scan Sequence Liabilities\" dialog should be unchecked", () => shouldBe(page, el("\"Free Cysteine\" checkbox in \"Scan Sequence Liabilities\" dialog"), "unchecked"));
      await session.step(51, "And \"Highlight in cell renderer\" checkbox in \"Scan Sequence Liabilities\" dialog should be checked", () => shouldBe(page, el("\"Highlight in cell renderer\" checkbox in \"Scan Sequence Liabilities\" dialog"), "checked"));
      await session.step(52, "And \"Create annotation column\" checkbox in \"Scan Sequence Liabilities\" dialog should be checked", () => shouldBe(page, el("\"Create annotation column\" checkbox in \"Scan Sequence Liabilities\" dialog"), "checked"));
      await session.step(53, "And \"Create summary count column\" checkbox in \"Scan Sequence Liabilities\" dialog should be unchecked", () => shouldBe(page, el("\"Create summary count column\" checkbox in \"Scan Sequence Liabilities\" dialog"), "unchecked"));
    });
    await run.scenario("The default rules write the per-row annotation column and hit real motifs", async () => {
      await session.step(56, "When user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
      await session.step(57, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(58, "And a new column \"~AntibodyHC_annotations\" should have been added", () => newColumnNamed(page, "~AntibodyHC_annotations"));
      await session.step(59, "And the table should not have a column \"AntibodyHC_liability_count\"", () => hasNoColumn(page, "AntibodyHC_liability_count"));
      await session.step(60, "And \"AntibodyHC\" column should carry at least 1 annotation", () => carriesAtLeast(page, "AntibodyHC", 1));
      await session.step(61, "And every liability hit on \"AntibodyHC\" column should match its motif at the position it reports", () => hitsMatch(page, "AntibodyHC"));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Manage Annotations lists one row per annotation and drops them one by one or all", async () => {
      await session.step(65, "When user picks \"Bio > Annotate > Manage Annotations...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Manage Annotations..."));
      await session.step(66, "Then \"Manage Annotations\" dialog should be visible", () => shouldBe(page, el("\"Manage Annotations\" dialog"), "visible"));
      await session.step(67, "And annotations list in \"Manage Annotations\" dialog should have as many items as \"AntibodyHC\" column has annotations", () => listMatchesAnnotations(page, "Manage Annotations", "AntibodyHC"));
      await session.step(68, "When user clicks on \"delete annotation\" icon in first item in annotations list", () => clickOn(page, el("\"delete annotation\" icon in first item in annotations list")));
      await session.step(69, "Then \"AntibodyHC\" column should carry one annotation fewer than before", () => oneFewer(page, "AntibodyHC"));
      await session.step(70, "And annotations list in \"Manage Annotations\" dialog should have as many items as \"AntibodyHC\" column has annotations", () => listMatchesAnnotations(page, "Manage Annotations", "AntibodyHC"));
      await session.step(71, "When user clicks on \"Clear All\" button in \"Manage Annotations\" dialog", () => clickOn(page, el("\"Clear All\" button in \"Manage Annotations\" dialog")));
      await session.step(72, "Then \"Manage Annotations\" dialog should contain text \"No annotations on this column.\"", () => shouldContainText(page, el("\"Manage Annotations\" dialog"), "No annotations on this column."));
      await session.step(73, "And \"AntibodyHC\" column should carry no annotations", () => carriesNone(page, "AntibodyHC"));
      await session.step(74, "When user clicks on CANCEL button in \"Manage Annotations\" dialog", () => clickOn(page, el("CANCEL button in \"Manage Annotations\" dialog")));
      await session.step(75, "Then \"Manage Annotations\" dialog should be hidden", () => shouldBe(page, el("\"Manage Annotations\" dialog"), "hidden"));
    });
    await run.scenario("Only the oxidation rules, reported as a count column", async () => {
      await session.step(78, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
      await session.step(79, "And user unchecks \"Deamidation (NG)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Deamidation (NG)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(80, "And user unchecks \"Deamidation (NS)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Deamidation (NS)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(81, "And user unchecks \"Deamidation (NA)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Deamidation (NA)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(82, "And user unchecks \"Deamidation (ND)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Deamidation (ND)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(83, "And user unchecks \"Deamidation (NT)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Deamidation (NT)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(84, "And user unchecks \"Isomerization (DG)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Isomerization (DG)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(85, "And user unchecks \"Isomerization (DS)\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Isomerization (DS)\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(86, "And user unchecks \"N-glycosylation\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"N-glycosylation\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(87, "And user unchecks \"Create annotation column\" checkbox in \"Scan Sequence Liabilities\" dialog", () => uncheck(page, el("\"Create annotation column\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(88, "And user checks \"Create summary count column\" checkbox in \"Scan Sequence Liabilities\" dialog", () => check(page, el("\"Create summary count column\" checkbox in \"Scan Sequence Liabilities\" dialog")));
      await session.step(89, "And user clicks on OK button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("OK button in \"Scan Sequence Liabilities\" dialog")));
      await session.step(90, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(91, "And a new column \"AntibodyHC_liability_count\" should have been added", () => newColumnNamed(page, "AntibodyHC_liability_count"));
      await session.step(92, "And \"AntibodyHC_liability_count\" column should have type \"int\"", () => columnType(page, "AntibodyHC_liability_count", "int"));
      await session.step(93, "And \"AntibodyHC_liability_count\" column should have no missing values", () => columnComplete(page, "AntibodyHC_liability_count"));
      await session.step(94, "And the total of \"AntibodyHC_liability_count\" column should be fewer than the liability hits found before", () => countBelowHits(page, "AntibodyHC_liability_count"));
      await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
