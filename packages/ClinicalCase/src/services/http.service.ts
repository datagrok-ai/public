import * as grok from "datagrok-api/grok";
import { CLINICAL_TRIAL_GOV_FIELDS } from "../constants/constants";

export const CLIN_TRIAL_GOV_SEARCH = 'https://clinicaltrials.gov/ct2/results?term=';

export class HttpService {

    clinTrialsGovUrl = 'https://clinicaltrials.gov/api/v2/';
    coverage = 'FullMatch';
    searchArea = 'OrgStudyId';
    format = 'json';

  constructor() {
  }

  async getStudyData(studyId: string, studyFields: string[]): Promise<{[key: string]: string} | null> {
    const url = `${this.clinTrialsGovUrl}studies?query.cond=COVERAGE%5B${this.coverage}%5DAREA%5B${this.searchArea}%5D${studyId}&fields=${studyFields.join(',')}&format=${this.format}`;
    const response = await grok.dapi.fetchProxy(url);
    const jsonResponse = await response.json();
    if (jsonResponse['studies'].length) {
      const response = jsonResponse['studies'][0]['protocolSection'];
      const res = {};
      res[CLINICAL_TRIAL_GOV_FIELDS.NCTId] = response['identificationModule']['nctId'];
      res[CLINICAL_TRIAL_GOV_FIELDS.Condition] = response['conditionsModule']['conditions'].join(',');
      res[CLINICAL_TRIAL_GOV_FIELDS.BriefTitle] = response['identificationModule']['briefTitle'];
      res[CLINICAL_TRIAL_GOV_FIELDS.OrgFullName] = response['identificationModule']['organization']['fullName'];
      res[CLINICAL_TRIAL_GOV_FIELDS.Phase] = response['designModule']['phases'].join(',');
      res[CLINICAL_TRIAL_GOV_FIELDS.StartDate] = response['statusModule']['startDateStruct']['date'];
      res[CLINICAL_TRIAL_GOV_FIELDS.CompletionDate] = response['statusModule']['completionDateStruct']['date'];
      res[CLINICAL_TRIAL_GOV_FIELDS.BriefSummary] = response['descriptionModule']['briefSummary'];
      return res;
    }
    return null;
  }

}