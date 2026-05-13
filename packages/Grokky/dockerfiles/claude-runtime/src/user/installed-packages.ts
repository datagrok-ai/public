import type {PackageInfo} from '../sync/packages';

const cache = new Map<string, Map<string, string>>();

export function getInstalledPackages(userId: string): Map<string, string> | undefined {
  return cache.get(userId);
}

export function isPackageInstalled(userId: string, name: string): boolean {
  return cache.get(userId)?.has(name) ?? false;
}

export function setInstalledPackages(userId: string, packages: PackageInfo[]): void {
  cache.set(userId, new Map(packages.map((p) => [p.name, p.updatedOn ?? ''])));
  console.log(`install-set: stored ${packages.length} package(s) for user ${userId}: [${packages.map((p) => p.name).join(', ')}]`);
}
