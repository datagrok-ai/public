export const CLIN_TRIAL_GOV_SEARCH = 'https://clinicaltrials.gov/ct2/results?term=';

export class HttpService {

    clinTrialsGovUrl = 'https://clinicaltrials.gov/api/query/';
    coverage = 'FullMatch';
    searchArea = 'OrgStudyId';
    format = 'JSON';

  constructor() {
  }

  async getStudyData(studyId: string, studyFields: string[]) {
    const response = await fetch(`${this.clinTrialsGovUrl}study_fields?expr=COVERAGE[${this.coverage}]AREA[${this.searchArea}]${studyId}&fields=${studyFields.join(',')}&fmt=${this.format}`);
    const jsonResponse = await response.json();
    if (jsonResponse['StudyFieldsResponse']['StudyFields']) {
      delete jsonResponse['StudyFieldsResponse']['StudyFields'][0]['Rank'];
      return jsonResponse['StudyFieldsResponse']['StudyFields'][0];
    }
    return null;
  }

}