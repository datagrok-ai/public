// Shared types for BenchlingLink package

export interface UserSummary {
  handle: string;
  id: string;
  name: string;
}

export interface Organization {
  handle: string;
  id: string;
  name: string;
}

export interface ArchiveRecord {
  reason?: string;
} 