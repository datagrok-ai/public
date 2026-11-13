/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/catalog.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../bindings/demo.js';
import {shouldBe, shouldHaveValue, toggle, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {el, enter, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Component catalog", () => {
  const session = feature(test, "features/demo/catalog.feature", import.meta.url);
  test("Component catalog", {tag: ["@demo", "@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 38, page);
    await session.step(10, "Given user opens the \"Catalog\" demo page", () => openDemoPage(page, "Catalog"));
    enter(page, "U2 Demo");
    await run.scenario("text input in every state", async () => {
      await session.step(13, "Then \"text-input-default\" text input should be enabled", () => shouldBe(page, el("\"text-input-default\" text input"), "enabled"));
      await session.step(14, "And \"text-input-disabled\" text input should be disabled", () => shouldBe(page, el("\"text-input-disabled\" text input"), "disabled"));
      await session.step(15, "And \"text-input-readOnly\" text input should be read-only", () => shouldBe(page, el("\"text-input-readOnly\" text input"), "read-only"));
      await session.step(16, "And \"text-input-invalid\" text input should be invalid", () => shouldBe(page, el("\"text-input-invalid\" text input"), "invalid"));
      await session.step(17, "When user types \"catalog\" into \"text-input-default\" text input", () => typeInto(page, "catalog", el("\"text-input-default\" text input")));
      await session.step(18, "Then \"text-input-default\" text input should have value \"catalog\"", () => shouldHaveValue(page, el("\"text-input-default\" text input"), "catalog"));
    });
    await run.scenario("text area in every state", async () => {
      await session.step(21, "Then \"text-area-default\" text area should be enabled", () => shouldBe(page, el("\"text-area-default\" text area"), "enabled"));
      await session.step(22, "And \"text-area-disabled\" text area should be disabled", () => shouldBe(page, el("\"text-area-disabled\" text area"), "disabled"));
      await session.step(23, "And \"text-area-readOnly\" text area should be read-only", () => shouldBe(page, el("\"text-area-readOnly\" text area"), "read-only"));
      await session.step(24, "And \"text-area-invalid\" text area should be invalid", () => shouldBe(page, el("\"text-area-invalid\" text area"), "invalid"));
      await session.step(25, "When user types \"catalog\" into \"text-area-default\" text area", () => typeInto(page, "catalog", el("\"text-area-default\" text area")));
      await session.step(26, "Then \"text-area-default\" text area should have value \"catalog\"", () => shouldHaveValue(page, el("\"text-area-default\" text area"), "catalog"));
    });
    await run.scenario("checkbox in every state", async () => {
      await session.step(29, "Then \"bool-input-default\" checkbox should be enabled", () => shouldBe(page, el("\"bool-input-default\" checkbox"), "enabled"));
      await session.step(30, "And \"bool-input-disabled\" checkbox should be disabled", () => shouldBe(page, el("\"bool-input-disabled\" checkbox"), "disabled"));
      await session.step(31, "And \"bool-input-readOnly\" checkbox should be read-only", () => shouldBe(page, el("\"bool-input-readOnly\" checkbox"), "read-only"));
      await session.step(32, "When user toggles \"bool-input-default\" checkbox", () => toggle(page, el("\"bool-input-default\" checkbox")));
      await session.step(33, "Then \"bool-input-default\" checkbox should be unchecked", () => shouldBe(page, el("\"bool-input-default\" checkbox"), "unchecked"));
    });
    await run.scenario("number input (u2-number-input) in every state", async () => {
      await session.step(36, "Then \"number-input-default\" number input should be enabled", () => shouldBe(page, el("\"number-input-default\" number input"), "enabled"));
      await session.step(37, "And \"number-input-disabled\" number input should be disabled", () => shouldBe(page, el("\"number-input-disabled\" number input"), "disabled"));
      await session.step(38, "And \"number-input-readOnly\" number input should be read-only", () => shouldBe(page, el("\"number-input-readOnly\" number input"), "read-only"));
      await session.step(39, "And \"number-input-invalid\" number input should be invalid", () => shouldBe(page, el("\"number-input-invalid\" number input"), "invalid"));
      await session.step(40, "When user types \"42\" into \"number-input-default\" number input", () => typeInto(page, "42", el("\"number-input-default\" number input")));
      await session.step(41, "Then \"number-input-default\" number input should have value \"42\"", () => shouldHaveValue(page, el("\"number-input-default\" number input"), "42"));
    });
    await run.scenario("choice input in every state", async () => {
      await session.step(44, "Then \"choice-input-default\" choice input should be enabled", () => shouldBe(page, el("\"choice-input-default\" choice input"), "enabled"));
      await session.step(45, "And \"choice-input-disabled\" choice input should be disabled", () => shouldBe(page, el("\"choice-input-disabled\" choice input"), "disabled"));
      await session.step(46, "And \"choice-input-readOnly\" choice input should be read-only", () => shouldBe(page, el("\"choice-input-readOnly\" choice input"), "read-only"));
      await session.step(47, "And \"choice-input-invalid\" choice input should be invalid", () => shouldBe(page, el("\"choice-input-invalid\" choice input"), "invalid"));
      await session.step(48, "And \"choice-input-empty\" choice input should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"choice-input-empty\" choice input"), 0));
    });
    await run.scenario("multi choice input in every state", async () => {
      await session.step(51, "Then \"multi-choice-input-default\" multi choice input should be enabled", () => shouldBe(page, el("\"multi-choice-input-default\" multi choice input"), "enabled"));
      await session.step(52, "And \"multi-choice-input-disabled\" multi choice input should be disabled", () => shouldBe(page, el("\"multi-choice-input-disabled\" multi choice input"), "disabled"));
      await session.step(53, "And \"multi-choice-input-readOnly\" multi choice input should be read-only", () => shouldBe(page, el("\"multi-choice-input-readOnly\" multi choice input"), "read-only"));
      await session.step(54, "And \"multi-choice-input-invalid\" multi choice input should be invalid", () => shouldBe(page, el("\"multi-choice-input-invalid\" multi choice input"), "invalid"));
      await session.step(55, "And \"multi-choice-input-empty\" multi choice input should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"multi-choice-input-empty\" multi choice input"), 0));
    });
    await run.scenario("button in every state", async () => {
      await session.step(58, "Then \"button-default\" button should be enabled", () => shouldBe(page, el("\"button-default\" button"), "enabled"));
      await session.step(59, "And \"button-disabled\" button should be disabled", () => shouldBe(page, el("\"button-disabled\" button"), "disabled"));
    });
    await run.scenario("horizontal stack in every state", async () => {
      await session.step(62, "Then \"form-buttons-default\" horizontal stack should be enabled", () => shouldBe(page, el("\"form-buttons-default\" horizontal stack"), "enabled"));
    });
    await run.scenario("splitter in every state", async () => {
      await session.step(65, "Then \"splitter-default\" splitter should be enabled", () => shouldBe(page, el("\"splitter-default\" splitter"), "enabled"));
    });
    await run.scenario("wizard in every state", async () => {
      await session.step(68, "Then \"wizard-default\" wizard should be enabled", () => shouldBe(page, el("\"wizard-default\" wizard"), "enabled"));
    });
    await run.scenario("accordion in every state", async () => {
      await session.step(71, "Then \"accordion-default\" accordion should be enabled", () => shouldBe(page, el("\"accordion-default\" accordion"), "enabled"));
    });
    await run.scenario("section in every state", async () => {
      await session.step(74, "Then \"section-default\" section should be enabled", () => shouldBe(page, el("\"section-default\" section"), "enabled"));
    });
    await run.scenario("tabs in every state", async () => {
      await session.step(77, "Then \"tabs-default\" tabs should be enabled", () => shouldBe(page, el("\"tabs-default\" tabs"), "enabled"));
    });
    await run.scenario("form in every state", async () => {
      await session.step(80, "Then \"form-default\" form should be enabled", () => shouldBe(page, el("\"form-default\" form"), "enabled"));
    });
    await run.scenario("property grid in every state", async () => {
      await session.step(83, "Then \"property-grid-default\" property grid should be enabled", () => shouldBe(page, el("\"property-grid-default\" property grid"), "enabled"));
    });
    await run.scenario("icon in every state", async () => {
      await session.step(86, "Then \"icon-default\" icon should be enabled", () => shouldBe(page, el("\"icon-default\" icon"), "enabled"));
    });
    await run.scenario("breadcrumbs in every state", async () => {
      await session.step(89, "Then \"breadcrumbs-default\" breadcrumbs should be enabled", () => shouldBe(page, el("\"breadcrumbs-default\" breadcrumbs"), "enabled"));
      await session.step(90, "And \"breadcrumbs-empty\" breadcrumbs should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"breadcrumbs-empty\" breadcrumbs"), 0));
    });
    await run.scenario("slider in every state", async () => {
      await session.step(93, "Then \"slider-input-default\" slider should be enabled", () => shouldBe(page, el("\"slider-input-default\" slider"), "enabled"));
      await session.step(94, "And \"slider-input-disabled\" slider should be disabled", () => shouldBe(page, el("\"slider-input-disabled\" slider"), "disabled"));
      await session.step(95, "And \"slider-input-readOnly\" slider should be read-only", () => shouldBe(page, el("\"slider-input-readOnly\" slider"), "read-only"));
      await session.step(96, "And \"slider-input-invalid\" slider should be invalid", () => shouldBe(page, el("\"slider-input-invalid\" slider"), "invalid"));
    });
    await run.scenario("range slider in every state", async () => {
      await session.step(99, "Then \"range-slider-default\" range slider should be enabled", () => shouldBe(page, el("\"range-slider-default\" range slider"), "enabled"));
    });
    await run.scenario("radio input in every state", async () => {
      await session.step(102, "Then \"radio-input-default\" radio input should be enabled", () => shouldBe(page, el("\"radio-input-default\" radio input"), "enabled"));
      await session.step(103, "And \"radio-input-disabled\" radio input should be disabled", () => shouldBe(page, el("\"radio-input-disabled\" radio input"), "disabled"));
      await session.step(104, "And \"radio-input-readOnly\" radio input should be read-only", () => shouldBe(page, el("\"radio-input-readOnly\" radio input"), "read-only"));
      await session.step(105, "And \"radio-input-invalid\" radio input should be invalid", () => shouldBe(page, el("\"radio-input-invalid\" radio input"), "invalid"));
      await session.step(106, "And \"radio-input-empty\" radio input should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"radio-input-empty\" radio input"), 0));
    });
    await run.scenario("color input in every state", async () => {
      await session.step(109, "Then \"color-input-default\" color input should be enabled", () => shouldBe(page, el("\"color-input-default\" color input"), "enabled"));
      await session.step(110, "And \"color-input-disabled\" color input should be disabled", () => shouldBe(page, el("\"color-input-disabled\" color input"), "disabled"));
      await session.step(111, "And \"color-input-readOnly\" color input should be read-only", () => shouldBe(page, el("\"color-input-readOnly\" color input"), "read-only"));
      await session.step(112, "When user types \"#ff0000\" into \"color-input-default\" color input", () => typeInto(page, "#ff0000", el("\"color-input-default\" color input")));
      await session.step(113, "Then \"color-input-default\" color input should have value \"#ff0000\"", () => shouldHaveValue(page, el("\"color-input-default\" color input"), "#ff0000"));
    });
    await run.scenario("icon input in every state", async () => {
      await session.step(116, "Then \"icon-input-default\" icon input should be enabled", () => shouldBe(page, el("\"icon-input-default\" icon input"), "enabled"));
      await session.step(117, "And \"icon-input-disabled\" icon input should be disabled", () => shouldBe(page, el("\"icon-input-disabled\" icon input"), "disabled"));
      await session.step(118, "And \"icon-input-readOnly\" icon input should be read-only", () => shouldBe(page, el("\"icon-input-readOnly\" icon input"), "read-only"));
      await session.step(119, "And \"icon-input-invalid\" icon input should be invalid", () => shouldBe(page, el("\"icon-input-invalid\" icon input"), "invalid"));
    });
    await run.scenario("list input in every state", async () => {
      await session.step(122, "Then \"list-input-default\" list input should be enabled", () => shouldBe(page, el("\"list-input-default\" list input"), "enabled"));
      await session.step(123, "And \"list-input-disabled\" list input should be disabled", () => shouldBe(page, el("\"list-input-disabled\" list input"), "disabled"));
      await session.step(124, "And \"list-input-readOnly\" list input should be read-only", () => shouldBe(page, el("\"list-input-readOnly\" list input"), "read-only"));
      await session.step(125, "And \"list-input-invalid\" list input should be invalid", () => shouldBe(page, el("\"list-input-invalid\" list input"), "invalid"));
    });
    await run.scenario("map input in every state", async () => {
      await session.step(128, "Then \"map-input-default\" map input should be enabled", () => shouldBe(page, el("\"map-input-default\" map input"), "enabled"));
      await session.step(129, "And \"map-input-disabled\" map input should be disabled", () => shouldBe(page, el("\"map-input-disabled\" map input"), "disabled"));
      await session.step(130, "And \"map-input-readOnly\" map input should be read-only", () => shouldBe(page, el("\"map-input-readOnly\" map input"), "read-only"));
      await session.step(131, "And \"map-input-invalid\" map input should be invalid", () => shouldBe(page, el("\"map-input-invalid\" map input"), "invalid"));
    });
    await run.scenario("records input in every state", async () => {
      await session.step(134, "Then \"records-input-default\" records input should be enabled", () => shouldBe(page, el("\"records-input-default\" records input"), "enabled"));
      await session.step(135, "And \"records-input-disabled\" records input should be disabled", () => shouldBe(page, el("\"records-input-disabled\" records input"), "disabled"));
      await session.step(136, "And \"records-input-readOnly\" records input should be read-only", () => shouldBe(page, el("\"records-input-readOnly\" records input"), "read-only"));
      await session.step(137, "And \"records-input-invalid\" records input should be invalid", () => shouldBe(page, el("\"records-input-invalid\" records input"), "invalid"));
    });
    await run.scenario("checklist in every state", async () => {
      await session.step(140, "Then \"checklist-input-default\" checklist should be enabled", () => shouldBe(page, el("\"checklist-input-default\" checklist"), "enabled"));
      await session.step(141, "And \"checklist-input-disabled\" checklist should be disabled", () => shouldBe(page, el("\"checklist-input-disabled\" checklist"), "disabled"));
      await session.step(142, "And \"checklist-input-readOnly\" checklist should be read-only", () => shouldBe(page, el("\"checklist-input-readOnly\" checklist"), "read-only"));
      await session.step(143, "And \"checklist-input-invalid\" checklist should be invalid", () => shouldBe(page, el("\"checklist-input-invalid\" checklist"), "invalid"));
      await session.step(144, "And \"checklist-input-empty\" checklist should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"checklist-input-empty\" checklist"), 0));
    });
    await run.scenario("chips input in every state", async () => {
      await session.step(147, "Then \"chips-input-default\" chips input should be enabled", () => shouldBe(page, el("\"chips-input-default\" chips input"), "enabled"));
      await session.step(148, "And \"chips-input-disabled\" chips input should be disabled", () => shouldBe(page, el("\"chips-input-disabled\" chips input"), "disabled"));
      await session.step(149, "And \"chips-input-readOnly\" chips input should be read-only", () => shouldBe(page, el("\"chips-input-readOnly\" chips input"), "read-only"));
      await session.step(150, "And \"chips-input-invalid\" chips input should be invalid", () => shouldBe(page, el("\"chips-input-invalid\" chips input"), "invalid"));
      await session.step(151, "And \"chips-input-empty\" chips input should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"chips-input-empty\" chips input"), 0));
    });
    await run.scenario("access grid in every state", async () => {
      await session.step(154, "Then \"access-grid-default\" access grid should be enabled", () => shouldBe(page, el("\"access-grid-default\" access grid"), "enabled"));
      await session.step(155, "And \"access-grid-disabled\" access grid should be disabled", () => shouldBe(page, el("\"access-grid-disabled\" access grid"), "disabled"));
      await session.step(156, "And \"access-grid-readOnly\" access grid should be read-only", () => shouldBe(page, el("\"access-grid-readOnly\" access grid"), "read-only"));
      await session.step(157, "And \"access-grid-invalid\" access grid should be invalid", () => shouldBe(page, el("\"access-grid-invalid\" access grid"), "invalid"));
    });
    await run.scenario("number input (u2-qnum-input) in every state", async () => {
      await session.step(160, "Then \"qnum-input-default\" number input should be enabled", () => shouldBe(page, el("\"qnum-input-default\" number input"), "enabled"));
      await session.step(161, "And \"qnum-input-disabled\" number input should be disabled", () => shouldBe(page, el("\"qnum-input-disabled\" number input"), "disabled"));
      await session.step(162, "And \"qnum-input-readOnly\" number input should be read-only", () => shouldBe(page, el("\"qnum-input-readOnly\" number input"), "read-only"));
      await session.step(163, "And \"qnum-input-invalid\" number input should be invalid", () => shouldBe(page, el("\"qnum-input-invalid\" number input"), "invalid"));
      await session.step(164, "When user types \"42\" into \"qnum-input-default\" number input", () => typeInto(page, "42", el("\"qnum-input-default\" number input")));
      await session.step(165, "Then \"qnum-input-default\" number input should have value \"42\"", () => shouldHaveValue(page, el("\"qnum-input-default\" number input"), "42"));
    });
    await run.scenario("font input in every state", async () => {
      await session.step(168, "Then \"font-input-default\" font input should be enabled", () => shouldBe(page, el("\"font-input-default\" font input"), "enabled"));
      await session.step(169, "And \"font-input-disabled\" font input should be disabled", () => shouldBe(page, el("\"font-input-disabled\" font input"), "disabled"));
      await session.step(170, "And \"font-input-readOnly\" font input should be read-only", () => shouldBe(page, el("\"font-input-readOnly\" font input"), "read-only"));
    });
    await run.scenario("image input in every state", async () => {
      await session.step(173, "Then \"image-input-default\" image input should be enabled", () => shouldBe(page, el("\"image-input-default\" image input"), "enabled"));
      await session.step(174, "And \"image-input-disabled\" image input should be disabled", () => shouldBe(page, el("\"image-input-disabled\" image input"), "disabled"));
      await session.step(175, "And \"image-input-readOnly\" image input should be read-only", () => shouldBe(page, el("\"image-input-readOnly\" image input"), "read-only"));
      await session.step(176, "And \"image-input-invalid\" image input should be invalid", () => shouldBe(page, el("\"image-input-invalid\" image input"), "invalid"));
    });
    await run.scenario("card in every state", async () => {
      await session.step(179, "Then \"card-default\" card should be enabled", () => shouldBe(page, el("\"card-default\" card"), "enabled"));
    });
    await run.scenario("stat card in every state", async () => {
      await session.step(182, "Then \"stat-card-default\" stat card should be enabled", () => shouldBe(page, el("\"stat-card-default\" stat card"), "enabled"));
    });
    await run.scenario("app header in every state", async () => {
      await session.step(185, "Then \"app-header-default\" app header should be enabled", () => shouldBe(page, el("\"app-header-default\" app header"), "enabled"));
      await session.step(186, "And \"app-header-empty\" app header should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"app-header-empty\" app header"), 0));
    });
    await run.scenario("table in every state", async () => {
      await session.step(189, "Then \"data-table-default\" table should be enabled", () => shouldBe(page, el("\"data-table-default\" table"), "enabled"));
      await session.step(190, "And \"data-table-empty\" table should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"data-table-empty\" table"), 0));
    });
    await run.scenario("catalog table in every state", async () => {
      await session.step(193, "Then \"catalog-table-default\" catalog table should be enabled", () => shouldBe(page, el("\"catalog-table-default\" catalog table"), "enabled"));
      await session.step(194, "And \"catalog-table-empty\" catalog table should have 0 items in total", () => shouldHaveItemsInTotal(page, el("\"catalog-table-empty\" catalog table"), 0));
    });
    await run.scenario("filter builder in every state", async () => {
      await session.step(197, "Then \"filter-builder-default\" filter builder should be enabled", () => shouldBe(page, el("\"filter-builder-default\" filter builder"), "enabled"));
      await session.step(198, "And \"filter-builder-disabled\" filter builder should be disabled", () => shouldBe(page, el("\"filter-builder-disabled\" filter builder"), "disabled"));
      await session.step(199, "And \"filter-builder-readOnly\" filter builder should be read-only", () => shouldBe(page, el("\"filter-builder-readOnly\" filter builder"), "read-only"));
      await session.step(200, "And \"filter-builder-invalid\" filter builder should be invalid", () => shouldBe(page, el("\"filter-builder-invalid\" filter builder"), "invalid"));
    });
    await run.scenario("filter query in every state", async () => {
      await session.step(203, "Then \"filter-query-input-default\" filter query should be enabled", () => shouldBe(page, el("\"filter-query-input-default\" filter query"), "enabled"));
      await session.step(204, "And \"filter-query-input-disabled\" filter query should be disabled", () => shouldBe(page, el("\"filter-query-input-disabled\" filter query"), "disabled"));
      await session.step(205, "And \"filter-query-input-readOnly\" filter query should be read-only", () => shouldBe(page, el("\"filter-query-input-readOnly\" filter query"), "read-only"));
      await session.step(206, "And \"filter-query-input-invalid\" filter query should be invalid", () => shouldBe(page, el("\"filter-query-input-invalid\" filter query"), "invalid"));
    });
    run.finish();
  });
});
