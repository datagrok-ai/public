import { study } from "../clinical-study";

export class ValidationHelper {

    domainsAndColumnsToCheck: {};
    missingColumnsInReqDomains;
    missingColumnsInOptDomains;
    missingDomains;

    constructor(domainsAndColumnsToCheck) {
        this.domainsAndColumnsToCheck = domainsAndColumnsToCheck;
        this.missingColumnsInReqDomains = {};
        this.missingColumnsInOptDomains = {};
        this.missingDomains = [];
    }

     private checkMissingDomains(requiredDomainsAndCols: any) {
        let reqDomains = requiredDomainsAndCols['req_domains'] ? Object.keys(requiredDomainsAndCols['req_domains']) : [];
        let optDomains = requiredDomainsAndCols['opt_domains'] ? Object.keys(requiredDomainsAndCols['opt_domains']) : [];
        let missingReqDomains = reqDomains.filter(it => study.domains[it] === null);
        let missingOptDomains = optDomains.some(it => study.domains[it] !== null) ? [] : optDomains;
        let presentOptDomains = optDomains.filter(it => study.domains[it] !== null);
        this.missingDomains = missingReqDomains.concat(missingOptDomains);
        let requiredColumns = {};
        reqDomains.forEach(domain => {
          requiredColumns[domain] = requiredDomainsAndCols['req_domains'][domain];
        });
        optDomains.forEach(domain => {
          requiredColumns[domain] = requiredDomainsAndCols['opt_domains'][domain];
        });
        if (!this.missingDomains.length) {
            reqDomains.forEach(it => {
                this.missingColumnsInReqDomains[it] = this.checkMissingColumns(it, requiredColumns);
            });
            presentOptDomains.forEach(it => {
                this.missingColumnsInOptDomains[it] = this.checkMissingColumns(it, requiredColumns);
            });
            
        } else {
            reqDomains.forEach(it => {
                this.missingColumnsInReqDomains[it] = this.checkMissingColumns(it, requiredColumns);
            });
            optDomains.forEach(it => {
                this.missingColumnsInOptDomains[it] = this.checkMissingColumns(it, requiredColumns);
            });
        }
      } 
      

      private checkMissingColumns(domain: string, requiredDomainsAndCols: any) {
          const domainColumns = study.domains[domain] ? study.domains[domain].columns.names() : [];
          const reqCols = requiredDomainsAndCols[domain]['req'] ?? [];
          const optCols = requiredDomainsAndCols[domain]['opt'] ?? []; //at least one of optional columns should exist in domain
          const missingReqColumns = reqCols.filter(it => !domainColumns.includes(it));
          const missingOptColumns = optCols.some(it => study.domains[it] !== null) ? [] : optCols.filter(it => !domainColumns.includes(it));
          const missingColumns = missingReqColumns.concat(missingOptColumns);
          return missingColumns;
      }


      validate(){
          if (!this.domainsAndColumnsToCheck) {
              return true;
          }
          this.checkMissingDomains(this.domainsAndColumnsToCheck);
          let missingReq = Object.keys(this.missingColumnsInReqDomains).length && Object.keys(this.missingColumnsInReqDomains).some(it => this.missingColumnsInReqDomains[it].length);
          let missingOpt = Object.keys(this.missingColumnsInOptDomains).length && Object.keys(this.missingColumnsInOptDomains).every(it => this.missingColumnsInOptDomains[it].length);
          return !this.missingDomains.length && !missingReq && !missingOpt;
      }

}