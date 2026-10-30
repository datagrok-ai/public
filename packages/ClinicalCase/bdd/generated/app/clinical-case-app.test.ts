/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/app/clinical-case-app.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, currentViewShows, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fileOnServer} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {errorBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Clinical Case opens from Browse > Apps with its demo studies and the Import study view", () => {
  const session = feature(test, "features/app/clinical-case-app.feature", import.meta.url);
  test("Clinical Case opens from Browse > Apps with its demo studies and the Import study view", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await run.scenario("The start page lists the demo studies with their subject counts", async () => {
      await session.step(22, "When user clicks on Apps---Clinical-Case tree node inside browse tree", () => clickOn(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(23, "Then the \"Clinical Case\" view should be current", () => viewIsCurrent(page, "Clinical Case"));
      await session.step(24, "And the current view should show the text \"Visualize and explore your SDTM data\"", () => currentViewShows(page, "Visualize and explore your SDTM data"));
      await session.step(25, "And the current view should show the text \"Study Standard Protocol Total subjects Start date End date CDISCPILOT01 SDTM 306 nida-ctn-0001 SDTM 138\"", () => currentViewShows(page, "Study Standard Protocol Total subjects Start date End date CDISCPILOT01 SDTM 306 nida-ctn-0001 SDTM 138"));
      await session.step(26, "And \"Import study...\" button should be visible", () => shouldBe(page, el("\"Import study...\" button"), "visible"));
      await session.step(28, "Given a file \"System:AppData/ClinicalCase/SDTM/CDISCPILOT01/studyConfig.json\" with text '{\"standard\":\"SDTM\",\"name\":\"CDISCPILOT01\",\"totalSubjects\":306}' is on the server", () => fileOnServer(page, "System:AppData/ClinicalCase/SDTM/CDISCPILOT01/studyConfig.json", "{\"standard\":\"SDTM\",\"name\":\"CDISCPILOT01\",\"totalSubjects\":306}"));
      await session.step(29, "And a file \"System:AppData/ClinicalCase/SDTM/nida-ctn-0001/studyConfig.json\" with text '{\"standard\":\"SDTM\",\"name\":\"nida-ctn-0001\",\"totalSubjects\":138,\"other\":{\"TEST PRODUCT IS BEING ADDED ON TO EXISTING TREATMENT\":\"Y\",\"AGE GROUP\":\"CHILDREN AND ADULTS\",\"PLANNED MAXIMUM AGE OF SUBJECTS\":\"NO UPPER LIMIT\",\"PLANNED MINIMUM AGE OF SUBJECTS\":\"15 YEARS\",\"COMPARATIVE TREATMENT NAME\":\"CLONIDINE\",\"TYPE OF CONTROL\":\"ACTIVE\",\"TRIAL DESIGN\":\"SINGLE GROUP PARALLEL\",\"DIAGNOSIS GROUP\":\"ADDICTED SUBJECTS\",\"TEST PRODUCT TOTAL DAILY DOSE\":\"ESCALATING FROM 4 MG BUP/1 MG NX TO 16 MG BUP TO 4 MG NX\",\"TEST PRODUCT DOSE UNITS\":\"MG\",\"DRUG DICTIONARY\":\"NONE\",\"TRIAL INDICATIONS\":\"OPIATE DEPENDENCE\",\"TRIAL INDICATION TYPE\":\"TREATMENT\",\"LENGTH OF TRIAL\":\"P6M\",\"TRIAL BLINDING SCHEMA\":\"OPEN LABEL\",\"TRIAL PRIMARY OBJECTIVE\":\"COMPARE THE RELATIVE CLINICAL UTILITY OF BUP/NX TO CLONIDINE IN A SHORT-TERM (13 DAY) DETOXIFICATION INITIATED IN AN INPATIENT SETTING\",\"TRIAL SECONDARY OBJECTIVE\":\"DETERMINE THE EFFECT OF TREATMENT ON WITHDRAWAL SYMPTOMS, TREATMENT RETENTION, SAFETY, PATIENT SATISFACTION, AND EXTENT OF ANCILLARY MEDICATIONS DISPENSED\",\"TRIAL PHASE\":\"4\",\"PLANNED NUMBER OF SUBJECTS\":\"360\",\"TRIAL IS RANDOMIZED\":\"Y\",\"RANDOMIZATION RATIO\":\"2:1 TEST TO COMPARATOR\",\"TEST PRODUCT ROUTE OF ADMINISTRATION\":\"SUBLINGUAL\",\"SEX OF PARTICIPANTS\":\"BOTH\",\"SPONSORING ORGANIZATION\":\"NIDA\",\"TRIAL TITLE\":\"BUPRENORPHINE/NALOXONE VERSUS CLONIDINE FOR INPATIENT OPIATE DETOXIFICATION\",\"REPORTED NAME OF TEST PRODUCT\":\"BUPRENORPHINE/CLONIDINE\",\"TYPE OF TRIAL\":\"EFFICACY\"}}' is on the server", () => fileOnServer(page, "System:AppData/ClinicalCase/SDTM/nida-ctn-0001/studyConfig.json", "{\"standard\":\"SDTM\",\"name\":\"nida-ctn-0001\",\"totalSubjects\":138,\"other\":{\"TEST PRODUCT IS BEING ADDED ON TO EXISTING TREATMENT\":\"Y\",\"AGE GROUP\":\"CHILDREN AND ADULTS\",\"PLANNED MAXIMUM AGE OF SUBJECTS\":\"NO UPPER LIMIT\",\"PLANNED MINIMUM AGE OF SUBJECTS\":\"15 YEARS\",\"COMPARATIVE TREATMENT NAME\":\"CLONIDINE\",\"TYPE OF CONTROL\":\"ACTIVE\",\"TRIAL DESIGN\":\"SINGLE GROUP PARALLEL\",\"DIAGNOSIS GROUP\":\"ADDICTED SUBJECTS\",\"TEST PRODUCT TOTAL DAILY DOSE\":\"ESCALATING FROM 4 MG BUP/1 MG NX TO 16 MG BUP TO 4 MG NX\",\"TEST PRODUCT DOSE UNITS\":\"MG\",\"DRUG DICTIONARY\":\"NONE\",\"TRIAL INDICATIONS\":\"OPIATE DEPENDENCE\",\"TRIAL INDICATION TYPE\":\"TREATMENT\",\"LENGTH OF TRIAL\":\"P6M\",\"TRIAL BLINDING SCHEMA\":\"OPEN LABEL\",\"TRIAL PRIMARY OBJECTIVE\":\"COMPARE THE RELATIVE CLINICAL UTILITY OF BUP/NX TO CLONIDINE IN A SHORT-TERM (13 DAY) DETOXIFICATION INITIATED IN AN INPATIENT SETTING\",\"TRIAL SECONDARY OBJECTIVE\":\"DETERMINE THE EFFECT OF TREATMENT ON WITHDRAWAL SYMPTOMS, TREATMENT RETENTION, SAFETY, PATIENT SATISFACTION, AND EXTENT OF ANCILLARY MEDICATIONS DISPENSED\",\"TRIAL PHASE\":\"4\",\"PLANNED NUMBER OF SUBJECTS\":\"360\",\"TRIAL IS RANDOMIZED\":\"Y\",\"RANDOMIZATION RATIO\":\"2:1 TEST TO COMPARATOR\",\"TEST PRODUCT ROUTE OF ADMINISTRATION\":\"SUBLINGUAL\",\"SEX OF PARTICIPANTS\":\"BOTH\",\"SPONSORING ORGANIZATION\":\"NIDA\",\"TRIAL TITLE\":\"BUPRENORPHINE/NALOXONE VERSUS CLONIDINE FOR INPATIENT OPIATE DETOXIFICATION\",\"REPORTED NAME OF TEST PRODUCT\":\"BUPRENORPHINE/CLONIDINE\",\"TYPE OF TRIAL\":\"EFFICACY\"}}"));
      await session.step(30, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Clinical Case node lists the studies and Import study", async () => {
      await session.step(33, "When user expands Apps---Clinical-Case tree node inside browse tree", () => expand(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(34, "Then Apps---Clinical-Case---CDISCPILOT01 tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---CDISCPILOT01 tree node inside browse tree"), "visible"));
      await session.step(35, "And Apps---Clinical-Case---nida-ctn-0001 tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---nida-ctn-0001 tree node inside browse tree"), "visible"));
      await session.step(36, "And Apps---Clinical-Case---Import-study tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---Import-study tree node inside browse tree"), "visible"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Import study... opens the import view, which refuses files without a study definition", async () => {
      await session.step(40, "When user clicks on Apps---Clinical-Case tree node inside browse tree", () => clickOn(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(41, "Then the \"Clinical Case\" view should be current", () => viewIsCurrent(page, "Clinical Case"));
      await session.step(42, "When user clicks on \"Import study...\" button", () => clickOn(page, el("\"Import study...\" button")));
      await session.step(43, "Then the \"Import Study - Clinical Case\" view should be current", () => viewIsCurrent(page, "Import Study - Clinical Case"));
      await session.step(44, "And the current view should show the text \"Import study files\"", () => currentViewShows(page, "Import study files"));
      await session.step(45, "And \"Import\" button should be disabled", () => shouldBe(page, el("\"Import\" button"), "disabled"));
      await session.step(46, "When user uploads \"../files/SDTM/CDISCPILOT01/dm.csv\" through \"Open files\" icon", () => uploadThrough(page, "../files/SDTM/CDISCPILOT01/dm.csv", el("\"Open files\" icon")));
      await session.step(47, "Then the current view should show the text \"Error: Invalid or missing define.xml/study.json\"", () => currentViewShows(page, "Error: Invalid or missing define.xml/study.json"));
      await session.step(48, "And the current view should show the text \"dm.csv\"", () => currentViewShows(page, "dm.csv"));
      await session.step(49, "And an error balloon containing \"Invalid or missing define.xml/study.json\" should have been shown", () => errorBalloonText(page, "Invalid or missing define.xml/study.json"));
      await session.step(50, "And \"Import\" button should be disabled", () => shouldBe(page, el("\"Import\" button"), "disabled"));
      await session.step(51, "When user closes the current view", () => closeCurrentView(page));
      await session.step(52, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
