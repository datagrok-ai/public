let _chemLocked = false;

export function chemLock() {
  if (_chemLocked) {
    throw('RdKit Service usage locked');
  }
  _chemLocked = true;
}

export function chemUnlock() {
  _chemLocked = false;
}
