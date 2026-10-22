export {DomainSource} from './core/domain-source.js';
export type {DomainSourceOptions, DomainSourceState, DomainSchema} from './core/domain-source.js';
export {DomainBackendError, emsBackends} from './core/domain-backend.js';
export type {DomainBackend, DomainTableLike, DomainTableInfoLike, DomainFrameLike, DomainQueryLike,
  DomainReadScope, DomainSupportLike, DomainTransactionOpLike, DomainTransactionResultLike,
  AuditEntryLike, DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchValidationLike,
  DomainDeletedMode, DomainProbeLike} from './core/domain-backend.js';
export {MemoryDomainBackend} from './core/memory-domain.js';
export type {MemorySchemaJson, MemoryTableJson, MemoryColumnJson, MemoryDomainOptions}
  from './core/memory-domain.js';
export type {EditState, EditSaved, RowState} from './core/edit-state.js';
export {SharedSession, confirmDiscard} from './core/session.js';
export type {DomainSession} from './core/session.js';
export {emsSourcesPlugin, registerDomainSources} from './core/registrations.js';
export {DomainQueries} from './core/queries.js';

export const VERSION = '0.1.0';
