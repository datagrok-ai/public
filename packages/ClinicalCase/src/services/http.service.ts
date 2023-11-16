import * as grok from "datagrok-api/grok";

export const CLIN_TRIAL_GOV_SEARCH = 'https://clinicaltrials.gov/ct2/results?term=';

export class HttpService {

    clinTrialsGovUrl = 'https://clinicaltrials.gov/api/query/';
    coverage = 'FullMatch';
    searchArea = 'OrgStudyId';
    format = 'JSON';

  constructor() {
  }

  async getStudyData(studyId: string, studyFields: string[]) {
    const url = `${this.clinTrialsGovUrl}study_fields?expr=COVERAGE%5B${this.coverage}%5DAREA%5B${this.searchArea}%5D${studyId}&fields=${studyFields.join(',')}&fmt=${this.format}`;
    const response = await grok.dapi.fetchProxy(url);
    const jsonResponse = await response.json();
    if (jsonResponse['StudyFieldsResponse']['StudyFields']) {
      delete jsonResponse['StudyFieldsResponse']['StudyFields'][0]['Rank'];
      return jsonResponse['StudyFieldsResponse']['StudyFields'][0];
    }
    return null;
  }

}