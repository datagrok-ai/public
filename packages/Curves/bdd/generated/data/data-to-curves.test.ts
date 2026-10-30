/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/data-to-curves.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {topMenuListsEnabled} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {call} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {autostartsCompleted, openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Data to Curves turns a table of wells into a table of fitted curves, one per compound, assay, target and run", () => {
  const session = feature(test, "features/data/data-to-curves.feature", import.meta.url);
  test("Data to Curves turns a table of wells into a table of fitted curves, one per compound, assay, target and run", {tag: ["@journey", "@known-failure:45"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(28, "And user opens a table \"wells\" with:", () => openTableOf(page, "wells", [["compound_id","assay_name","target_entity","batch_id","run_id","concentration_uM","inhibition"],["C1","A","T","B1","R1","0.1","5"],["C1","A","T","B1","R1","1","20"],["C1","A","T","B1","R1","10","60"],["C1","A","T","B1","R1","100","95"],["C1","A","T","B1","R2","0.1","4"],["C1","A","T","B1","R2","1","25"],["C1","A","T","B1","R2","10","55"],["C1","A","T","B1","R2","100","90"],["C2","A","T","B2","R1","0.1","2"],["C2","A","T","B2","R1","1","10"],["C2","A","T","B2","R1","10","30"],["C2","A","T","B2","R1","100","50"]]), [["compound_id","assay_name","target_entity","batch_id","run_id","concentration_uM","inhibition"],["C1","A","T","B1","R1","0.1","5"],["C1","A","T","B1","R1","1","20"],["C1","A","T","B1","R1","10","60"],["C1","A","T","B1","R1","100","95"],["C1","A","T","B1","R2","0.1","4"],["C1","A","T","B1","R2","1","25"],["C1","A","T","B1","R2","10","55"],["C1","A","T","B1","R2","100","90"],["C2","A","T","B2","R1","0.1","2"],["C2","A","T","B2","R1","1","10"],["C2","A","T","B2","R1","10","30"],["C2","A","T","B2","R1","100","50"]]);
    await run.scenario("The Data menu offers Data to Curves", async () => {
      await session.step(45, "Then the top menu should list \"Data > Curves > Data to Curves\" as enabled", () => topMenuListsEnabled(page, "Data > Curves > Data to Curves"));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true, line: 45});
    await run.scenario("The dialog prefills the roles from the column names, and OK opens the Fitted Curves view", async () => {
      await session.step(49, "When user calls \"Curves:dataToCurvesTopMenu\" function", () => call(page, "Curves:dataToCurvesTopMenu"));
      await session.step(50, "Then \"Data to Curves\" dialog should be visible", () => shouldBe(page, el("\"Data to Curves\" dialog"), "visible"));
      await session.step(51, "And \"Well Table\" input in \"Data to Curves\" dialog should contain text \"wells\"", () => shouldContainText(page, el("\"Well Table\" input in \"Data to Curves\" dialog"), "wells"));
      await session.step(52, "And Concentration input in \"Data to Curves\" dialog should contain text \"concentration_uM\"", () => shouldContainText(page, el("Concentration input in \"Data to Curves\" dialog"), "concentration_uM"));
      await session.step(53, "And Readout input in \"Data to Curves\" dialog should contain text \"inhibition\"", () => shouldContainText(page, el("Readout input in \"Data to Curves\" dialog"), "inhibition"));
      await session.step(54, "And \"Batch ID\" input in \"Data to Curves\" dialog should contain text \"batch_id\"", () => shouldContainText(page, el("\"Batch ID\" input in \"Data to Curves\" dialog"), "batch_id"));
      await session.step(55, "And \"Run ID\" input in \"Data to Curves\" dialog should contain text \"run_id\"", () => shouldContainText(page, el("\"Run ID\" input in \"Data to Curves\" dialog"), "run_id"));
      await session.step(56, "And \"Compound ID\" input in \"Data to Curves\" dialog should contain text \"compound_id\"", () => shouldContainText(page, el("\"Compound ID\" input in \"Data to Curves\" dialog"), "compound_id"));
      await session.step(57, "And \"Assay Name\" input in \"Data to Curves\" dialog should contain text \"assay_name\"", () => shouldContainText(page, el("\"Assay Name\" input in \"Data to Curves\" dialog"), "assay_name"));
      await session.step(58, "And \"Target Entity\" input in \"Data to Curves\" dialog should contain text \"target_entity\"", () => shouldContainText(page, el("\"Target Entity\" input in \"Data to Curves\" dialog"), "target_entity"));
      await session.step(59, "When user clicks on OK button in \"Data to Curves\" dialog", () => clickOn(page, el("OK button in \"Data to Curves\" dialog")));
      await session.step(60, "Then \"Data to Curves\" dialog should be hidden", () => shouldBe(page, el("\"Data to Curves\" dialog"), "hidden"));
      await session.step(61, "And the open table views should be exactly \"wells, Fitted Curves\"", () => openTableViewsExactly(page, "wells, Fitted Curves"));
      await session.step(62, "Given user switches to the \"Fitted Curves\" table view", () => switchTableView(page, "Fitted Curves"));
      await session.step(63, "Then the table should have 3 rows", () => rowCount(page, 3));
      await session.step(64, "And \"Fitted Curve\" column should have semantic type \"fit\"", () => columnSemType(page, "Fitted Curve", "fit"));
      await session.step(65, "And the value of \"Compound|Assay|Target\" column in row 1 should be \"C1 | A | T\"", () => valueInRow(page, "Compound|Assay|Target", 1, "C1 | A | T"));
      await session.step(66, "And the value of \"Run ID\" column in row 1 should be \"R1\"", () => valueInRow(page, "Run ID", 1, "R1"));
      await session.step(67, "And the value of \"Max Percent Inhibition\" column in row 1 should be \"95\"", () => valueInRow(page, "Max Percent Inhibition", 1, "95"));
      await session.step(68, "And the value of \"Compound|Assay|Target\" column in row 2 should be \"C1 | A | T\"", () => valueInRow(page, "Compound|Assay|Target", 2, "C1 | A | T"));
      await session.step(69, "And the value of \"Run ID\" column in row 2 should be \"R2\"", () => valueInRow(page, "Run ID", 2, "R2"));
      await session.step(70, "And the value of \"Max Percent Inhibition\" column in row 2 should be \"90\"", () => valueInRow(page, "Max Percent Inhibition", 2, "90"));
      await session.step(71, "And the value of \"Compound|Assay|Target\" column in row 3 should be \"C2 | A | T\"", () => valueInRow(page, "Compound|Assay|Target", 3, "C2 | A | T"));
      await session.step(72, "And the value of \"Max Percent Inhibition\" column in row 3 should be \"50\"", () => valueInRow(page, "Max Percent Inhibition", 3, "50"));
      await session.step(73, "And trellis plot viewer should be visible", () => shouldBe(page, el("trellis plot viewer"), "visible"));
      await session.step(74, "And the open tableview should have 1 MultiCurveViewer viewer", () => viewerCount(page, 1, "MultiCurveViewer"));
      await session.step(75, "And table \"wells\" should have 12 rows", () => tableRows(page, "wells", 12));
      await session.step(76, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
