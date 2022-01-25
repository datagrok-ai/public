import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {delay, expect} from "@datagrok-libraries/utils/src/test";
import {_package} from "../package-test";
import { StudySummaryView } from "../views/study-summary-view";
import { ClinicalDomains, study } from "../clinical-study";
import { createBaselineEndpointDataframe, createHysLawDataframe, createLabValuesByVisitDataframe, createSurvivalData, cumulativeEnrollemntByDay, dynamicComparedToBaseline, labDynamicComparedToMinMax, labDynamicRelatedToRef } from "../data-preparation/data-preparation";
import { AE_START_DATE, LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, SUBJECT_ID, SUBJ_REF_ENDT, SUBJ_REF_STDT, TREATMENT_ARM, VISIT_DAY, VISIT_NAME } from "../columns-constants";
import { dataframeContentToRow } from "../data-preparation/utils";
import { createValidationDataFrame } from "../validation/validation-utils";
import { vaidateDMDomain } from "../validation/services/validation-service";
import { ALT, BILIRUBIN } from "../constants";
import { StudyVisit } from "../model/study-visit";
import { PatientVisit } from "../model/patient-visit";

export const allViews = [
  'Summary', 
  'Timelines', 
  'Patient Profile', 
  'Adverse Events', 
  'Laboratory', 
  'Survival Analysis', 
  'Distributions', 
  'Correlations', 
  'Time Profile', 
  'Tree map', 
  'Medical History', 
  'Visits'];

  export const allTableViews = ['AE Browser'];

export async function requireText(name: string): Promise<string> {
    return await _package.files.readAsText(name);
  }

export async function _testOpenApp() {
  await createTableView('dm.csv');
  await grok.functions.call("Clinicalcase:clinicalCaseApp");
  await delay(1000);
  expect(grok.shell.v.name === 'Summary', true);
}

export async function _testAllViewsCreated() { 
  expect(allViews.every(item => grok.shell.view(item) !== undefined), true);
}

export async function _testAllTableViewsCreated() { 
  expect(allTableViews.every(item => grok.shell.view(item) !== undefined), true);
}

export async function _testCumulativeEnrollment() {
  const df = await readDataframe('dm.csv');
  const cumulativeEnrDf = cumulativeEnrollemntByDay(df, SUBJ_REF_STDT, SUBJECT_ID, 'CUMULATIVE_ENROLLMENT');
  expect(cumulativeEnrDf.rowCount, 5);
  expect(cumulativeEnrDf.get(SUBJ_REF_STDT, 0).toString(), '2012-08-04 23:00:00.000');
  expect(cumulativeEnrDf.get('CUMULATIVE_ENROLLMENT', 4), 6);
}

export async function _testValidation() {
  const df = await readDataframe('dm.csv');
  const validationDf = createValidationDataFrame();
  vaidateDMDomain(df, validationDf);
  expect(validationDf.rowCount, 4);
  testRowValues(validationDf, 0, {'Column':'RFSTDTC', 'Row number': 6, 'Value': '', 'Violated rule ID': 'SD0003'});
}

export async function _testHysLaw() {
  const lb = await readDataframe('lb.csv');
  const dm = await readDataframe('dm.csv');
  const hysLawDf = createHysLawDataframe(lb, dm, 'Alanine Aminotransferase', 'Aspartate Aminotransferase', 'Bilirubin');
  expect(hysLawDf.rowCount, 6);
  testRowValues(hysLawDf, 0, {[SUBJECT_ID]: '01-701-1015', [ALT]: 3.5, [BILIRUBIN]: 3});
}

export async function _testBaselineEndpoint() {
  const lb = await readDataframe('lb.csv');
  const dm = await readDataframe('dm.csv');
  const baselineEndpointDataframe = createBaselineEndpointDataframe(lb, dm, [TREATMENT_ARM], LAB_TEST, LAB_RES_N, 
    [LAB_LO_LIM_N, LAB_HI_LIM_N], 'Alkaline Phosphatase', 'SCREENING 1', 'WEEK 4', VISIT_NAME, 'BL', 'EP');
  expect(baselineEndpointDataframe.rowCount, 6);
  testRowValues(baselineEndpointDataframe, 0, {[SUBJECT_ID]: '01-701-1015', 'BL': 34, 'EP': 41});
}

export async function _testLabValuesByVisit() {
  const lb = await readDataframe('lb.csv');
  const dm = await readDataframe('dm.csv');
  const labByVisitDataframe = createLabValuesByVisitDataframe(lb, dm, 'Alkaline Phosphatase', 
    'Placebo', 'Values', VISIT_DAY);
  expect(labByVisitDataframe.rowCount, 40);
  testRowValues(labByVisitDataframe, 0, {[SUBJECT_ID]: '01-701-1015', 'Values': 34, [VISIT_DAY]: -7});
}

export async function _testSurvivalDataframe() {
  const ae = await readDataframe('ae.csv');
  const dm = await readDataframe('dm.csv');
  const survivalDataframe = createSurvivalData(dm, null, 'RETAIN IN STUDY', SUBJ_REF_ENDT, ['SEX']);
  const survivalDataframeDrugRel = createSurvivalData(dm, ae, 'DRUG RELATED AE', AE_START_DATE, ['SEX']);
  expect(survivalDataframe.rowCount, 6);
  expect(survivalDataframeDrugRel.rowCount, 6);
  testRowValues(survivalDataframe, 0, {[SUBJECT_ID]: '01-701-1015', 'time': 181, 'status': 1, 'SEX': 'F'});
  testRowValues(survivalDataframeDrugRel, 0, {[SUBJECT_ID]: '01-701-1015', 'time': 1, 'status': 1, 'SEX': 'F'});
}

export async function _testLabChangesComparedToBaseline(){
  const lb = await readDataframe('lb.csv');
  dynamicComparedToBaseline(lb, LAB_TEST, LAB_RES_N, 'SCREENING 1', VISIT_NAME, 'LAB_DYNAMIC_BL', false);
  expect(lb.rowCount, 1329);
  testRowValues(lb, 1, {[SUBJECT_ID]: '01-701-1015', 'LBTEST': 'Albumin', 'BL_LBSTRESN': 38, [LAB_RES_N]: 76, 'LAB_DYNAMIC_BL': 1});
}

export async function _testLabChangesBetweenMinMax(){
  const lb = await readDataframe('lb.csv');
  labDynamicComparedToMinMax(lb, 'LAB_DYNAMIC_MIN_MAX');
  expect(lb.rowCount, 1329);
  testRowValues(lb, 3, {[SUBJECT_ID]: '01-701-1015', 'LBTEST': 'Albumin', 'min(LBSTRESN)': 16.5, 'max(LBSTRESN)': 98, [LAB_RES_N]: 98, 'LAB_DYNAMIC_MIN_MAX': 1});
}

export async function _testLabChangesRelatedToRef(){
  const lb = await readDataframe('lb.csv');
  labDynamicRelatedToRef(lb, 'LAB_DYNAMIC_REF');
  expect(lb.rowCount, 1329);
  testRowValues(lb, 0, {[SUBJECT_ID]: '01-701-1015', 'LBTEST': 'Albumin', [LAB_LO_LIM_N]: 33, [LAB_HI_LIM_N]: 49, [LAB_RES_N]: 38, 'LAB_DYNAMIC_REF': -0.375});
  testRowValues(lb, 2, {[SUBJECT_ID]: '01-701-1015', 'LBTEST': 'Albumin', [LAB_LO_LIM_N]: 33, [LAB_HI_LIM_N]: 49, [LAB_RES_N]: 16.5, 'LAB_DYNAMIC_REF': -2});
  testRowValues(lb, 3, {[SUBJECT_ID]: '01-701-1015', 'LBTEST': 'Albumin', [LAB_LO_LIM_N]: 33, [LAB_HI_LIM_N]: 49, [LAB_RES_N]: 98, 'LAB_DYNAMIC_REF': 2});
}

export async function _testStudyVisit(){
  const domains = await createClinicalDomains(['sv', 'lb', 'ae']);
  const studyVisit = new StudyVisit(domains);
  studyVisit.updateStudyVisit(28, 'WEEK 4', 14);
  expect(studyVisit.totalPatients, 6);
  expect(studyVisit.minVisitDate, '02.09.2012');
  expect(studyVisit.maxVisitDate, '29.07.2014');
  expect(studyVisit.lbAtVisit.rowCount, 191);
  expect(studyVisit.aeSincePreviusVisit.rowCount, 4);
}

export async function _testPatientVisit(){
  const domains = await createClinicalDomains(['sv', 'lb', 'ae']);
  const patientVisit = new PatientVisit(domains);
  patientVisit.updateSubjectVisit('01-701-1028', 28, 'WEEK 4', 14);
  patientVisit.updateSubjectVisitDomains();
  console.log(dataframeContentToRow(patientVisit['lb']));
  console.log(dataframeContentToRow(patientVisit['ae']));
  expect(patientVisit['lb'].rowCount, 30);
  expect(patientVisit['ae'].rowCount, 1);
}

async function createClinicalDomains(domainNames: string[]){
  const domains = new ClinicalDomains();
  await Promise.all(domainNames.map(async (it) => {
    const df = await readDataframe(`${it}.csv`);
    domains[it] = df;
  }));
  return domains;
}

async function createTableView(tableName: string) {
  const df = await readDataframe(tableName);
  df.name = tableName.replace('.csv', '');
  grok.shell.addTableView(df);
}

async function readDataframe(tableName: string) {
  const file = await requireText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}

function testRowValues(df: DG.DataFrame, rowNumber: number, columnsToCheck: any) {
  Object.keys(columnsToCheck).forEach(key => expect(df.get(key, rowNumber), columnsToCheck[key]));
}
