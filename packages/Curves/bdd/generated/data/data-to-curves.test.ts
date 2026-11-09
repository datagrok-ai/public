/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/data-to-curves.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [core#1246]
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data to Curves turns a table of wells into a table of fitted curves, one per compound, assay, target and run", () => {
  const session = feature(test, "features/data/data-to-curves.feature", import.meta.url);
  test("The Data menu opens the dialog with the roles prefilled from the column names, and OK opens the Fitted Curves view", {tag: ["@realizes:core#1246"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(23, "And user opens a table \"wells\" with:", () => openTableOf(page, "wells", [["compound_id","assay_name","target_entity","batch_id","run_id","concentration_uM","inhibition"],["C1","A","T","B1","R1","0.1","5"],["C1","A","T","B1","R1","1","20"],["C1","A","T","B1","R1","10","60"],["C1","A","T","B1","R1","100","95"],["C1","A","T","B1","R2","0.1","4"],["C1","A","T","B1","R2","1","25"],["C1","A","T","B1","R2","10","55"],["C1","A","T","B1","R2","100","90"],["C2","A","T","B2","R1","0.1","2"],["C2","A","T","B2","R1","1","10"],["C2","A","T","B2","R1","10","30"],["C2","A","T","B2","R1","100","50"]]), [["compound_id","assay_name","target_entity","batch_id","run_id","concentration_uM","inhibition"],["C1","A","T","B1","R1","0.1","5"],["C1","A","T","B1","R1","1","20"],["C1","A","T","B1","R1","10","60"],["C1","A","T","B1","R1","100","95"],["C1","A","T","B1","R2","0.1","4"],["C1","A","T","B1","R2","1","25"],["C1","A","T","B1","R2","10","55"],["C1","A","T","B1","R2","100","90"],["C2","A","T","B2","R1","0.1","2"],["C2","A","T","B2","R1","1","10"],["C2","A","T","B2","R1","10","30"],["C2","A","T","B2","R1","100","50"]]);
    await session.step(39, "When user picks \"Data > Curves > Data to Curves\" from the top menu", () => pickFromTopMenu(page, "Data > Curves > Data to Curves"));
    await session.step(40, "Then \"Data to Curves\" dialog should be visible", () => shouldBe(page, el("\"Data to Curves\" dialog"), "visible"));
    await session.step(41, "And Concentration input in \"Data to Curves\" dialog should contain text \"concentration_uM\"", () => shouldContainText(page, el("Concentration input in \"Data to Curves\" dialog"), "concentration_uM"));
    await session.step(42, "And Readout input in \"Data to Curves\" dialog should contain text \"inhibition\"", () => shouldContainText(page, el("Readout input in \"Data to Curves\" dialog"), "inhibition"));
    await session.step(43, "When user clicks on OK button in \"Data to Curves\" dialog", () => clickOn(page, el("OK button in \"Data to Curves\" dialog")));
    await session.step(44, "Then the open table views should be exactly \"wells, Fitted Curves\"", () => openTableViewsExactly(page, "wells, Fitted Curves"));
    await session.step(45, "Given user switches to the \"Fitted Curves\" table view", () => switchTableView(page, "Fitted Curves"));
    await session.step(46, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(47, "And \"Fitted Curve\" column should have semantic type \"fit\"", () => columnSemType(page, "Fitted Curve", "fit"));
    await session.step(48, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
    await session.step(49, "And the open tableview should have 1 MultiCurveViewer viewer", () => viewerCount(page, 1, "MultiCurveViewer"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
  });
});
