/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/menu/top-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.top-menu.registration]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {monomerSketcherReady} from '../../bindings/bio-b.js';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldHaveText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, everyValueContains, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {commandCompleted, newColumnMatching, pickFromTopMenu, topMenuLists} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterIsExactlyContains, filterPassesAll, filterPassesFewer, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDatasetRowsAs, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, painted, propertyShouldBe, readingAtLeast} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Bio top menu", () => {
  const session = feature(test, "features/menu/top-menu.feature", import.meta.url);
  test("The Bio top menu", {tag: ["@journey", "@realizes:bio.top-menu.registration"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 25, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens antibodies dataset keeping the first 40 rows as \"antibodies\"", () => openDatasetRowsAs(page, ds("antibodies"), 40, "antibodies"));
    await session.step(19, "And the Bio package is initialized", () => bioInitialized(page));
    await run.scenario("Every group lists its commands", async () => {
      await session.step(22, "Then the top menu should list:", () => topMenuLists(page, [["Bio > Transform > Molecules to HELM..."],["Bio > Transform > To Atomic Level..."],["Bio > Transform > Convert Sequence Notation..."],["Bio > Transform > Split to Monomers..."],["Bio > Analyze > Activity Cliffs..."],["Bio > Analyze > Sequence Space..."],["Bio > Analyze > MSA..."],["Bio > Analyze > Compare sequences..."],["Bio > Analyze > Composition"],["Bio > Calculate > Extract Region..."],["Bio > Calculate > Identity..."],["Bio > Calculate > Similarity..."],["Bio > Annotate > Apply Numbering Scheme..."],["Bio > Annotate > Scan Liabilities..."],["Bio > Annotate > Manage Annotations..."],["Bio > Manage > Match with Monomer Library..."],["Bio > Manage > Monomer Libraries"],["Bio > Manage > Monomers"],["Bio > Search > Similarity Search"],["Bio > Search > Diversity Search"],["Bio > Search > Subsequence Search ..."]]), [["Bio > Transform > Molecules to HELM..."],["Bio > Transform > To Atomic Level..."],["Bio > Transform > Convert Sequence Notation..."],["Bio > Transform > Split to Monomers..."],["Bio > Analyze > Activity Cliffs..."],["Bio > Analyze > Sequence Space..."],["Bio > Analyze > MSA..."],["Bio > Analyze > Compare sequences..."],["Bio > Analyze > Composition"],["Bio > Calculate > Extract Region..."],["Bio > Calculate > Identity..."],["Bio > Calculate > Similarity..."],["Bio > Annotate > Apply Numbering Scheme..."],["Bio > Annotate > Scan Liabilities..."],["Bio > Annotate > Manage Annotations..."],["Bio > Manage > Match with Monomer Library..."],["Bio > Manage > Monomer Libraries"],["Bio > Manage > Monomers"],["Bio > Search > Similarity Search"],["Bio > Search > Diversity Search"],["Bio > Search > Subsequence Search ..."]]);
    });
    await run.scenario("Bio > Transform > Molecules to HELM... opens the Molecules to HELM dialog and cancels cleanly [group=Transform, leaf=Molecules to HELM..., dialog=Molecules to HELM]", async () => {
      await session.step(46, "When user picks \"Bio > Transform > Molecules to HELM...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Molecules to HELM..."));
      await session.step(47, "Then \"Molecules to HELM\" dialog should be visible", () => shouldBe(page, el("\"Molecules to HELM\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Molecules to HELM\" dialog should be visible", () => shouldBe(page, el("OK button in \"Molecules to HELM\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Molecules to HELM\" dialog", () => clickOn(page, el("CANCEL button in \"Molecules to HELM\" dialog")));
      await session.step(50, "Then \"Molecules to HELM\" dialog should be hidden", () => shouldBe(page, el("\"Molecules to HELM\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Transform > To Atomic Level... opens the To Atomic Level dialog and cancels cleanly [group=Transform, leaf=To Atomic Level..., dialog=To Atomic Level]", async () => {
      await session.step(46, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(47, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(48, "And OK button in \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("OK button in \"To Atomic Level\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"To Atomic Level\" dialog", () => clickOn(page, el("CANCEL button in \"To Atomic Level\" dialog")));
      await session.step(50, "Then \"To Atomic Level\" dialog should be hidden", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Transform > Convert Sequence Notation... opens the Convert Sequence Notation dialog and cancels cleanly [group=Transform, leaf=Convert Sequence Notation..., dialog=Convert Sequence Notation]", async () => {
      await session.step(46, "When user picks \"Bio > Transform > Convert Sequence Notation...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Convert Sequence Notation..."));
      await session.step(47, "Then \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Convert Sequence Notation\" dialog should be visible", () => shouldBe(page, el("OK button in \"Convert Sequence Notation\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Convert Sequence Notation\" dialog", () => clickOn(page, el("CANCEL button in \"Convert Sequence Notation\" dialog")));
      await session.step(50, "Then \"Convert Sequence Notation\" dialog should be hidden", () => shouldBe(page, el("\"Convert Sequence Notation\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Transform > Split to Monomers... opens the Split to Monomers dialog and cancels cleanly [group=Transform, leaf=Split to Monomers..., dialog=Split to Monomers]", async () => {
      await session.step(46, "When user picks \"Bio > Transform > Split to Monomers...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > Split to Monomers..."));
      await session.step(47, "Then \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Split to Monomers\" dialog should be visible", () => shouldBe(page, el("OK button in \"Split to Monomers\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Split to Monomers\" dialog", () => clickOn(page, el("CANCEL button in \"Split to Monomers\" dialog")));
      await session.step(50, "Then \"Split to Monomers\" dialog should be hidden", () => shouldBe(page, el("\"Split to Monomers\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Analyze > Activity Cliffs... opens the Sequence Activity Cliffs dialog and cancels cleanly [group=Analyze, leaf=Activity Cliffs..., dialog=Sequence Activity Cliffs]", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > Activity Cliffs...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Activity Cliffs..."));
      await session.step(47, "Then \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Sequence Activity Cliffs\" dialog should be visible", () => shouldBe(page, el("OK button in \"Sequence Activity Cliffs\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Sequence Activity Cliffs\" dialog", () => clickOn(page, el("CANCEL button in \"Sequence Activity Cliffs\" dialog")));
      await session.step(50, "Then \"Sequence Activity Cliffs\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Activity Cliffs\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Analyze > Sequence Space... opens the Sequence Space dialog and cancels cleanly [group=Analyze, leaf=Sequence Space..., dialog=Sequence Space]", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > Sequence Space...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Sequence Space..."));
      await session.step(47, "Then \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("\"Sequence Space\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Sequence Space\" dialog should be visible", () => shouldBe(page, el("OK button in \"Sequence Space\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Sequence Space\" dialog", () => clickOn(page, el("CANCEL button in \"Sequence Space\" dialog")));
      await session.step(50, "Then \"Sequence Space\" dialog should be hidden", () => shouldBe(page, el("\"Sequence Space\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Analyze > MSA... opens the MSA dialog and cancels cleanly [group=Analyze, leaf=MSA..., dialog=MSA]", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > MSA...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > MSA..."));
      await session.step(47, "Then \"MSA\" dialog should be visible", () => shouldBe(page, el("\"MSA\" dialog"), "visible"));
      await session.step(48, "And OK button in \"MSA\" dialog should be visible", () => shouldBe(page, el("OK button in \"MSA\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"MSA\" dialog", () => clickOn(page, el("CANCEL button in \"MSA\" dialog")));
      await session.step(50, "Then \"MSA\" dialog should be hidden", () => shouldBe(page, el("\"MSA\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Analyze > Compare sequences... opens the Compare Sequences dialog and cancels cleanly [group=Analyze, leaf=Compare sequences..., dialog=Compare Sequences]", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > Compare sequences...\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Compare sequences..."));
      await session.step(47, "Then \"Compare Sequences\" dialog should be visible", () => shouldBe(page, el("\"Compare Sequences\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Compare Sequences\" dialog should be visible", () => shouldBe(page, el("OK button in \"Compare Sequences\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Compare Sequences\" dialog", () => clickOn(page, el("CANCEL button in \"Compare Sequences\" dialog")));
      await session.step(50, "Then \"Compare Sequences\" dialog should be hidden", () => shouldBe(page, el("\"Compare Sequences\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Analyze > Composition opens the Composition Analysis dialog and cancels cleanly [group=Analyze, leaf=Composition, dialog=Composition Analysis]", async () => {
      await session.step(46, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
      await session.step(47, "Then \"Composition Analysis\" dialog should be visible", () => shouldBe(page, el("\"Composition Analysis\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Composition Analysis\" dialog should be visible", () => shouldBe(page, el("OK button in \"Composition Analysis\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Composition Analysis\" dialog", () => clickOn(page, el("CANCEL button in \"Composition Analysis\" dialog")));
      await session.step(50, "Then \"Composition Analysis\" dialog should be hidden", () => shouldBe(page, el("\"Composition Analysis\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Calculate > Extract Region... opens the Get Sequence Region dialog and cancels cleanly [group=Calculate, leaf=Extract Region..., dialog=Get Sequence Region]", async () => {
      await session.step(46, "When user picks \"Bio > Calculate > Extract Region...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Extract Region..."));
      await session.step(47, "Then \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Get Sequence Region\" dialog should be visible", () => shouldBe(page, el("OK button in \"Get Sequence Region\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Get Sequence Region\" dialog", () => clickOn(page, el("CANCEL button in \"Get Sequence Region\" dialog")));
      await session.step(50, "Then \"Get Sequence Region\" dialog should be hidden", () => shouldBe(page, el("\"Get Sequence Region\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Calculate > Identity... opens the Identity dialog and cancels cleanly [group=Calculate, leaf=Identity..., dialog=Identity]", async () => {
      await session.step(46, "When user picks \"Bio > Calculate > Identity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Identity..."));
      await session.step(47, "Then \"Identity\" dialog should be visible", () => shouldBe(page, el("\"Identity\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Identity\" dialog should be visible", () => shouldBe(page, el("OK button in \"Identity\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Identity\" dialog", () => clickOn(page, el("CANCEL button in \"Identity\" dialog")));
      await session.step(50, "Then \"Identity\" dialog should be hidden", () => shouldBe(page, el("\"Identity\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Calculate > Similarity... opens the Similarity dialog and cancels cleanly [group=Calculate, leaf=Similarity..., dialog=Similarity]", async () => {
      await session.step(46, "When user picks \"Bio > Calculate > Similarity...\" from the top menu", () => pickFromTopMenu(page, "Bio > Calculate > Similarity..."));
      await session.step(47, "Then \"Similarity\" dialog should be visible", () => shouldBe(page, el("\"Similarity\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Similarity\" dialog should be visible", () => shouldBe(page, el("OK button in \"Similarity\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Similarity\" dialog", () => clickOn(page, el("CANCEL button in \"Similarity\" dialog")));
      await session.step(50, "Then \"Similarity\" dialog should be hidden", () => shouldBe(page, el("\"Similarity\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Annotate > Apply Numbering Scheme... opens the Apply Antibody Numbering dialog and cancels cleanly [group=Annotate, leaf=Apply Numbering Scheme..., dialog=Apply Antibody Numbering]", async () => {
      await session.step(46, "When user picks \"Bio > Annotate > Apply Numbering Scheme...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Apply Numbering Scheme..."));
      await session.step(47, "Then \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Apply Antibody Numbering\" dialog should be visible", () => shouldBe(page, el("OK button in \"Apply Antibody Numbering\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Apply Antibody Numbering\" dialog", () => clickOn(page, el("CANCEL button in \"Apply Antibody Numbering\" dialog")));
      await session.step(50, "Then \"Apply Antibody Numbering\" dialog should be hidden", () => shouldBe(page, el("\"Apply Antibody Numbering\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Annotate > Scan Liabilities... opens the Scan Sequence Liabilities dialog and cancels cleanly [group=Annotate, leaf=Scan Liabilities..., dialog=Scan Sequence Liabilities]", async () => {
      await session.step(46, "When user picks \"Bio > Annotate > Scan Liabilities...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Scan Liabilities..."));
      await session.step(47, "Then \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Scan Sequence Liabilities\" dialog should be visible", () => shouldBe(page, el("OK button in \"Scan Sequence Liabilities\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Scan Sequence Liabilities\" dialog", () => clickOn(page, el("CANCEL button in \"Scan Sequence Liabilities\" dialog")));
      await session.step(50, "Then \"Scan Sequence Liabilities\" dialog should be hidden", () => shouldBe(page, el("\"Scan Sequence Liabilities\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Annotate > Manage Annotations... opens the Manage Annotations dialog and cancels cleanly [group=Annotate, leaf=Manage Annotations..., dialog=Manage Annotations]", async () => {
      await session.step(46, "When user picks \"Bio > Annotate > Manage Annotations...\" from the top menu", () => pickFromTopMenu(page, "Bio > Annotate > Manage Annotations..."));
      await session.step(47, "Then \"Manage Annotations\" dialog should be visible", () => shouldBe(page, el("\"Manage Annotations\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Manage Annotations\" dialog should be visible", () => shouldBe(page, el("OK button in \"Manage Annotations\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Manage Annotations\" dialog", () => clickOn(page, el("CANCEL button in \"Manage Annotations\" dialog")));
      await session.step(50, "Then \"Manage Annotations\" dialog should be hidden", () => shouldBe(page, el("\"Manage Annotations\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Manage > Match with Monomer Library... opens the Match with Monomer Library dialog and cancels cleanly [group=Manage, leaf=Match with Monomer Library..., dialog=Match with Monomer Library]", async () => {
      await session.step(46, "When user picks \"Bio > Manage > Match with Monomer Library...\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Match with Monomer Library..."));
      await session.step(47, "Then \"Match with Monomer Library\" dialog should be visible", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Match with Monomer Library\" dialog should be visible", () => shouldBe(page, el("OK button in \"Match with Monomer Library\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Match with Monomer Library\" dialog", () => clickOn(page, el("CANCEL button in \"Match with Monomer Library\" dialog")));
      await session.step(50, "Then \"Match with Monomer Library\" dialog should be hidden", () => shouldBe(page, el("\"Match with Monomer Library\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Search > Subsequence Search ... opens the Substructure Search dialog and cancels cleanly [group=Search, leaf=Subsequence Search ..., dialog=Substructure Search]", async () => {
      await session.step(46, "When user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
      await session.step(47, "Then \"Substructure Search\" dialog should be visible", () => shouldBe(page, el("\"Substructure Search\" dialog"), "visible"));
      await session.step(48, "And OK button in \"Substructure Search\" dialog should be visible", () => shouldBe(page, el("OK button in \"Substructure Search\" dialog"), "visible"));
      await session.step(49, "When user clicks on CANCEL button in \"Substructure Search\" dialog", () => clickOn(page, el("CANCEL button in \"Substructure Search\" dialog")));
      await session.step(50, "Then \"Substructure Search\" dialog should be hidden", () => shouldBe(page, el("\"Substructure Search\" dialog"), "hidden"));
      await session.step(51, "And dialog should be hidden", () => shouldBe(page, el("dialog"), "hidden"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Bio > Manage > Monomer Libraries opens the manager with the libraries listed", async () => {
      await session.step(74, "When user picks \"Bio > Manage > Monomer Libraries\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomer Libraries"));
      await session.step(75, "Then the \"Manage Monomer Libraries\" view should be current", () => viewIsCurrent(page, "Manage Monomer Libraries"));
      await session.step(76, "And \"HELMCoreLibrary.json\" checkbox should be visible", () => shouldBe(page, el("\"HELMCoreLibrary.json\" checkbox"), "visible"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
      await session.step(78, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Bio > Manage > Monomers opens the monomer table", async () => {
      await session.step(81, "When user picks \"Bio > Manage > Monomers\" from the top menu", () => pickFromTopMenu(page, "Bio > Manage > Monomers"));
      await session.step(82, "Then the \"Manage Monomers\" view should be current", () => viewIsCurrent(page, "Manage Monomers"));
      await session.step(83, "And the table should have a column \"Symbol\"", () => hasColumn(page, "Symbol"));
      await session.step(84, "And \"Symbol\" column should have no missing values", () => columnComplete(page, "Symbol"));
      await session.step(85, "And the monomer sketcher of the Manage Monomers view should be ready", () => monomerSketcherReady(page));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
      await session.step(87, "When user closes the current view", () => closeCurrentView(page));
    });
    await run.scenario("Bio > Search > Similarity Search docks a \"Sequence Similarity Search\" viewer that has computed [leaf=Similarity Search, viewer=Sequence Similarity Search, reading=neighbours]", async () => {
      await session.step(90, "When user picks \"Bio > Search > Similarity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Similarity Search"));
      await session.step(91, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(92, "And \"Sequence Similarity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Similarity Search\" viewer"), "visible"));
      await session.step(93, "And the \"neighbours\" reading of \"Sequence Similarity Search\" viewer should be at least 2", () => readingAtLeast(page, "neighbours", el("\"Sequence Similarity Search\" viewer"), 2));
    });
    await run.scenario("Bio > Search > Diversity Search docks a \"Sequence Diversity Search\" viewer that has computed [leaf=Diversity Search, viewer=Sequence Diversity Search, reading=subset size]", async () => {
      await session.step(90, "When user picks \"Bio > Search > Diversity Search\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Diversity Search"));
      await session.step(91, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(92, "And \"Sequence Diversity Search\" viewer should be visible", () => shouldBe(page, el("\"Sequence Diversity Search\" viewer"), "visible"));
      await session.step(93, "And the \"subset size\" reading of \"Sequence Diversity Search\" viewer should be at least 2", () => readingAtLeast(page, "subset size", el("\"Sequence Diversity Search\" viewer"), 2));
    });
    await run.scenario("Subsequence Search on two sequence columns asks which, and its OK filters by the query", async () => {
      await session.step(100, "When user switches to the \"antibodies\" table view", () => switchTableView(page, "antibodies"));
      await session.step(101, "And user picks \"Bio > Search > Subsequence Search ...\" from the top menu", () => pickFromTopMenu(page, "Bio > Search > Subsequence Search ..."));
      await session.step(102, "Then \"Substructure Search\" dialog should be visible", () => shouldBe(page, el("\"Substructure Search\" dialog"), "visible"));
      await session.step(103, "And editor of Column input in \"Substructure Search\" dialog should have text \"AntibodyHC\"", () => shouldHaveText(page, el("editor of Column input in \"Substructure Search\" dialog"), "AntibodyHC"));
      await session.step(104, "When user enters \"SCAASGFTINGT\" into Substructure input in \"Substructure Search\" dialog", () => enterInto(page, "SCAASGFTINGT", el("Substructure input in \"Substructure Search\" dialog")));
      await session.step(105, "And user clicks on OK button in \"Substructure Search\" dialog", () => clickOn(page, el("OK button in \"Substructure Search\" dialog")));
      await session.step(106, "Then \"Substructure Search\" dialog should be hidden", () => shouldBe(page, el("\"Substructure Search\" dialog"), "hidden"));
      await session.step(107, "And fewer than 40 rows should pass the filter", () => filterPassesFewer(page, 40));
      await session.step(108, "And the filter should pass exactly the rows where \"AntibodyHC\" contains \"SCAASGFTINGT\"", () => filterIsExactlyContains(page, "AntibodyHC", "SCAASGFTINGT"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
      await session.step(110, "When user resets the filter", () => resetFilter(page));
      await session.step(111, "Then all rows should pass the filter", () => filterPassesAll(page));
    });
    await run.scenario("Composition on the second sequence column binds its WebLogo there", async () => {
      await session.step(114, "When user picks \"Bio > Analyze > Composition\" from the top menu", () => pickFromTopMenu(page, "Bio > Analyze > Composition"));
      await session.step(115, "Then \"Composition Analysis\" dialog should be visible", () => shouldBe(page, el("\"Composition Analysis\" dialog"), "visible"));
      await session.step(116, "When user selects \"AntibodyLC\" in Column input in \"Composition Analysis\" dialog", () => selectIn(page, "AntibodyLC", el("Column input in \"Composition Analysis\" dialog")));
      await session.step(117, "And user clicks on OK button in \"Composition Analysis\" dialog", () => clickOn(page, el("OK button in \"Composition Analysis\" dialog")));
      await session.step(118, "Then WebLogo viewer should be visible", () => shouldBe(page, el("WebLogo viewer"), "visible"));
      await session.step(119, "And \"Sequence Column Name\" property of WebLogo viewer should be \"AntibodyLC\"", () => propertyShouldBe(page, "Sequence Column Name", el("WebLogo viewer"), "AntibodyLC"));
      await session.step(120, "And WebLogo viewer should be painted", () => painted(page, el("WebLogo viewer")));
      await session.step(121, "And no errors should have been logged", () => noErrors(page));
      await session.step(122, "When user clicks on close icon of WebLogo viewer", () => clickOn(page, el("close icon of WebLogo viewer")));
      await session.step(123, "Then WebLogo viewer should be absent", () => shouldBe(page, el("WebLogo viewer"), "absent"));
    });
    await run.scenario("To Atomic Level on the second column with Non-linear off still converts every row", async () => {
      await session.step(126, "When user picks \"Bio > Transform > To Atomic Level...\" from the top menu", () => pickFromTopMenu(page, "Bio > Transform > To Atomic Level..."));
      await session.step(127, "Then \"To Atomic Level\" dialog should be visible", () => shouldBe(page, el("\"To Atomic Level\" dialog"), "visible"));
      await session.step(128, "When user selects \"AntibodyLC\" in Sequence input in \"To Atomic Level\" dialog", () => selectIn(page, "AntibodyLC", el("Sequence input in \"To Atomic Level\" dialog")));
      await session.step(129, "And user unchecks \"Non-linear\" checkbox in \"To Atomic Level\" dialog", () => uncheck(page, el("\"Non-linear\" checkbox in \"To Atomic Level\" dialog")));
      await session.step(130, "And user clicks on OK button in \"To Atomic Level\" dialog", () => clickOn(page, el("OK button in \"To Atomic Level\" dialog")));
      await session.step(131, "Then the top menu command should have completed", () => commandCompleted(page));
      await session.step(132, "And a new column matching \"^molfile\\(AntibodyLC\\)\" should have been added", () => newColumnMatching(page, "^molfile\\(AntibodyLC\\)"));
      await session.step(133, "And \"molfile(AntibodyLC)\" column should have no missing values", () => columnComplete(page, "molfile(AntibodyLC)"));
      await session.step(134, "And every value of \"molfile(AntibodyLC)\" column should contain \"M  V30 BEGIN CTAB\"", () => everyValueContains(page, "molfile(AntibodyLC)", "M  V30 BEGIN CTAB"));
      await session.step(135, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
