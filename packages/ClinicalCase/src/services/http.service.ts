export const CLIN_TRIAL_GOV_SEARCH = 'https://clinicaltrials.gov/ct2/results?term=';
export const NCTID_FIELD = 'NCTId';

export class HttpService {

    clinTrialsGovUrl = 'https://clinicaltrials.gov/api/query/';
    studyFields = [NCTID_FIELD,'Condition','BriefTitle','OrgFullName','Phase','StartDate','CompletionDate','BriefSummary'];
    coverage = 'FullMatch';
    searchArea = 'OrgStudyId';
    format = 'JSON';

  constructor() {
  }

  async getStudyData(studyId: string) {
      const response = await fetch(`${this.clinTrialsGovUrl}study_fields?expr=COVERAGE[${this.coverage}]AREA[${this.searchArea}]${studyId}&fields=${this.studyFields.join(',')}&fmt=${this.format}`);
      const jsonResponse =  await response.json();
      delete jsonResponse['StudyFieldsResponse']['StudyFields'][0]['Rank'];
      return jsonResponse['StudyFieldsResponse']['StudyFields'][0];
  }

}