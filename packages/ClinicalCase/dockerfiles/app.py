import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import pyreadstat
from celery import Celery
from filelock import FileLock
from datagrok_celery_task import DatagrokTask, Settings
from datagrok_api import DatagrokClient


settings = Settings(log_level=logging.DEBUG)
app = Celery(settings.celery_name, broker=settings.broker_url)

CORE_BINARY = Path(os.environ.get("CORE_BINARY", "/app/core"))
DEFAULT_TIMEOUT = int(os.environ.get("CORE_TIMEOUT_SECONDS", "1800"))
DATASETS_ROOT = Path(os.environ.get("DATASETS_ROOT", "/app/datasets"))
LOCKS_ROOT = Path(os.environ.get("DATASETS_LOCKS_ROOT", "/tmp/datasets_locks"))
DATASETS_ROOT.mkdir(parents=True, exist_ok=True)
LOCKS_ROOT.mkdir(parents=True, exist_ok=True)


def _ensure_core_exists() -> None:
    if not CORE_BINARY.exists():
        raise RuntimeError(f"CDISC CORE binary is missing at {CORE_BINARY}.")


def _quoted(cmd: List[str]) -> str:
    return " ".join(shlex.quote(arg) for arg in cmd)


def _run_core(arguments: List[str], timeout: Optional[int] = None) -> subprocess.CompletedProcess[str]:
    _ensure_core_exists()
    cmd = [str(CORE_BINARY), *arguments]
    logging.info("Running command: %s", _quoted(cmd))
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout or DEFAULT_TIMEOUT,
    )
    # Log RSS so we can spot OOM pressure in retrospect — when the kernel OOM-kills CORE,
    # subprocess.run returns no info; this baseline at least confirms whether our own process
    # is approaching the container memory limit.
    try:
        with open("/proc/self/status", "r") as f:
            for line in f:
                if line.startswith("VmPeak:") or line.startswith("VmRSS:"):
                    logging.info("worker mem %s", line.strip())
    except OSError:
        pass
    logging.debug("Command finished with code %s", proc.returncode)
    return proc


def _append_flag(args: List[str], flag: str, value: Optional[str]) -> None:
    if value:
        args.extend([flag, value])


def _append_repeatable_flag(args: List[str], flag: str, values: Optional[List[str]]) -> None:
    if values:
        for val in values:
            if val:
                args.extend([flag, val])


def _build_validation_args(
    standard: str,
    version: Optional[str],
    data_path: str,
    output_format: str,
    options: Optional[Dict[str, object]] = None,
) -> List[str]:
    if not standard:
        raise RuntimeError("Standard is required.")
    if not data_path:
        raise RuntimeError("Data path is required.")

    args: List[str] = [
        "validate",
        "-s",
        standard,
    ]
    _append_flag(args, "--version", version)
    args.extend([
        "-d",
        data_path,
        "-of",
        output_format or "json",
    ])

    opts = options or {}
    _append_flag(args, "--cache-path", opts.get("cache_path"))
    _append_flag(args, "--custom-standard", opts.get("custom_standard"))
    _append_flag(args, "--custom-standard-version", opts.get("custom_standard_version"))
    _append_flag(args, "--custom-standard-file", opts.get("custom_standard_file"))
    _append_flag(args, "--datasets", opts.get("datasets"))
    _append_repeatable_flag(args, "-r", opts.get("rule_ids"))

    extra_args = opts.get("extra_args")
    if isinstance(extra_args, list):
        args.extend(str(arg) for arg in extra_args if arg)
    elif isinstance(extra_args, str):
        args.extend(shlex.split(extra_args))

    return args


def _csv_to_xpt(csv_path: Path, xpt_path: Path) -> None:
    """Convert a single CSV file to SAS XPT format suitable for CDISC CORE validation."""
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False, na_values=[""])
    # XPT v5 requires uppercase column names of <=8 chars — SDTM/SEND conform by spec,
    # but be defensive against truncation collisions: log if any occur after normalization.
    normalized = [re.sub(r"[^A-Z0-9_]", "_", c.upper())[:8] for c in df.columns]
    if len(set(normalized)) != len(normalized):
        raise RuntimeError(
            f"Cannot convert {csv_path.name}: column-name normalization produced duplicates "
            f"(original: {list(df.columns)}, normalized: {normalized})"
        )
    df.columns = normalized
    # Try to coerce numeric-looking columns back to numeric so XPT stores them as numeric, not character.
    for col in df.columns:
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.notna().sum() == df[col].notna().sum() and df[col].notna().any():
            df[col] = coerced
    # Table name in XPT is also limited to 8 uppercase chars.
    table_name = re.sub(r"[^A-Z0-9_]", "_", csv_path.stem.upper())[:8]
    # CDISC CORE expects SAS Transport v5 (the legacy XPT format used by SDTM/SEND submissions).
    pyreadstat.write_xport(df, str(xpt_path), table_name=table_name, file_format_version=5)


def _convert_csvs_to_xpt_dir(src_dir: Path) -> Optional[Path]:
    """If src_dir contains any .csv datasets, mirror its contents into a temp dir
    with each .csv converted to .xpt (and existing .xpt / define.xml passed through).
    Returns the temp dir path, or None if no conversion was needed.
    The caller is responsible for cleaning up the returned directory.
    """
    csvs = [p for p in src_dir.iterdir() if p.is_file() and p.suffix.lower() == ".csv"]
    if not csvs:
        return None

    converted_dir = Path(tempfile.mkdtemp(prefix="xpt_", dir=str(src_dir.parent)))
    logging.info(f"Converting CSV datasets in {src_dir} to XPT in {converted_dir}")

    # Pass through only files CORE understands: existing .xpt datasets and define.xml.
    # Everything else (Datagrok caches like .d42, .sync_meta.json[.backup.NN], studyConfig.json,
    # arbitrary .csv post-processing artifacts, etc.) is intentionally excluded.
    passthrough_exts = {".xpt"}
    passthrough_names = {"define.xml"}
    for item in src_dir.iterdir():
        if not item.is_file():
            continue
        if item.suffix.lower() in passthrough_exts or item.name.lower() in passthrough_names:
            shutil.copy2(item, converted_dir / item.name)

    for csv in csvs:
        xpt_path = converted_dir / f"{csv.stem}.xpt"
        # If a real .xpt with the same stem was already copied through, prefer it.
        if xpt_path.exists():
            logging.info(f"Skipping conversion of {csv.name}: {xpt_path.name} already present")
            continue
        try:
            _csv_to_xpt(csv, xpt_path)
            logging.info(f"Converted {csv.name} -> {xpt_path.name}")
        except Exception as e:
            logging.error(f"Failed to convert {csv.name}: {e}")
            raise

    return converted_dir


def _write_records_as_csv(records: List[Dict[str, object]], dest: Path) -> None:
    """Write a list of dict records as CSV (UTF-8, comma-separated, RFC 4180 quoting).
    Nested arrays / dicts become JSON-encoded strings inside the cell so JS can
    JSON.parse them on demand. Empty/missing input writes a header-only file with no rows
    (DG.DataFrame.fromCsv handles that gracefully)."""
    import csv as _csv
    # Union of keys across all records — handles records with different shapes.
    keys: List[str] = []
    seen = set()
    for r in records:
        if not isinstance(r, dict):
            continue
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(dest, "w", encoding="utf-8", newline="") as f:
        writer = _csv.writer(f, quoting=_csv.QUOTE_MINIMAL)
        writer.writerow(keys)
        for r in records:
            if not isinstance(r, dict):
                continue
            row = []
            for k in keys:
                v = r.get(k)
                if v is None:
                    row.append("")
                elif isinstance(v, (list, dict)):
                    row.append(json.dumps(v, separators=(",", ":")))
                else:
                    row.append(v)
            writer.writerow(row)


def _sync_directory_with(api: DatagrokClient, remote_path: str, dest: Path) -> None:
    logging.info("Syncing %s to %s", remote_path, dest)
    api.files.sync_dir("System:AppData", remote_path, str(dest), recursive=True)
    logging.info("Sync complete: %s", dest)


#name: Run CDISC CORE Validation
#input: string standard
#input: string version { optional: true }
#input: string data_path
#input: string output_format { optional: true }
#input: string options { optional: true }
#output: string result
@app.task(name="run_core_validate", bind=True, base=DatagrokTask)
def run_core_validate(
    self,
    standard: str,
    version: Optional[str],
    data_path: str,
    output_format: str = "json",
    options: Optional[str] = None,
    **kwargs,
) -> str:
    self.update_state(meta={"description": "Preparing CDISC CORE validation."})
    
    if not data_path:
        raise RuntimeError("data_path is required.")

    # Check if data_path is a remote path (starts with System:AppData/)
    token = kwargs.get("USER_API_KEY")
    if not token:
        raise RuntimeError("USER_API_KEY is required to sync remote data.")
    if not settings.api_url:
        raise RuntimeError("DATAGROK_API_URL is not configured for the container.")
    api = DatagrokClient(api_key=token, base_url=settings.api_url)
    local_dir = DATASETS_ROOT / data_path
    logging.info(f"local_dir {local_dir}")
    lock_file = LOCKS_ROOT / f"{data_path.replace('/', '_')}.lock"
    logging.info(f"lock_file {lock_file}")
    with FileLock(str(lock_file), timeout=240):
        _sync_directory_with(api, data_path, local_dir)
    local_data_path = str(local_dir)

    # Temporarily move .sync_meta.json file if it exists (created by Datagrok file sync)
    # Validation fails if folder contains JSON files that aren't study datasets
    sync_meta_file = Path(local_data_path) / ".sync_meta.json"
    sync_meta_backup: Optional[Path] = None
    if sync_meta_file.exists():
        try:
            sync_meta_backup = sync_meta_file.parent / f".sync_meta.json.backup.{os.getpid()}"
            sync_meta_file.rename(sync_meta_backup)
            logging.info(f"Temporarily moved .sync_meta.json to {sync_meta_backup}")
        except OSError as e:
            logging.warning(f"Failed to move .sync_meta.json: {e}")

    # Remove studyConfig.json file if it exists
    # Validation fails if folder contains JSON files that aren't study datasets
    study_config_file = Path(local_data_path) / "studyConfig.json"
    if study_config_file.exists():
        try:
            study_config_file.unlink()
            logging.info(f"Removed studyConfig.json from {local_data_path}")
        except OSError as e:
            logging.warning(f"Failed to remove studyConfig.json: {e}")

    converted_dir: Optional[Path] = None
    try:
        # Parse options JSON string if provided
        options_dict: Optional[Dict[str, object]] = None
        if options:
            try:
                options_dict = json.loads(options)
            except json.JSONDecodeError as e:
                raise RuntimeError(f"Invalid JSON in options parameter: {e}")

        # CDISC CORE only accepts SAS Transport (.xpt) datasets — convert any CSVs first.
        self.update_state(meta={"description": "Converting CSV datasets to XPT."})
        converted_dir = _convert_csvs_to_xpt_dir(Path(local_data_path))
        validation_path = str(converted_dir) if converted_dir else local_data_path

        logging.info(f"validation_path {validation_path}")
        try:
            dir_listing = sorted(p.name for p in Path(validation_path).iterdir() if p.is_file())
            logging.info(f"validation_path contents ({len(dir_listing)} files): {dir_listing}")
        except OSError as e:
            logging.warning(f"Could not list validation_path: {e}")

        args = _build_validation_args(standard, version, validation_path, output_format, options_dict)

        self.update_state(meta={"description": "Running validation inside container."})
        proc = _run_core(args)

        logging.info(f"core exit code: {proc.returncode}")
        logging.info(f"core stdout: {proc.stdout!r}")
        logging.info(f"core stderr: {proc.stderr!r}")

        if proc.returncode not in (0, 1):
            raise RuntimeError(
                f"core validate failed with exit code {proc.returncode}. stderr: {proc.stderr.strip()}"
            )

        # proc.stdout contains the path to the result JSON file in format: "Output: CORE-Report-2025-12-01T18-51-22"
        stdout_text = proc.stdout.strip()
        if not stdout_text:
            raise RuntimeError(
                f"Validation completed but no result file path in stdout. "
                f"exit_code={proc.returncode}, stderr: {proc.stderr.strip() or '<empty>'}"
            )

        # Extract filename after "Output: "
        if not stdout_text.startswith("Output: "):
            raise RuntimeError(f"Unexpected stdout format: {stdout_text}")
        
        result_file_name = stdout_text[len("Output: "):].strip()
        if not result_file_name:
            raise RuntimeError("No filename found after 'Output: ' in stdout")

        # CORE versions differ: older builds print the filename without an extension,
        # newer ones include it. Try the name as-is first, then with the format extension appended.
        ext = (output_format or "json").lstrip(".")
        candidates = [Path(f"/app/{result_file_name}"), Path(f"/app/{result_file_name}.{ext}")]
        result_file = next((p for p in candidates if p.exists()), None)
        if result_file is None:
            raise RuntimeError(
                f"Result file not found. Tried: {[str(p) for p in candidates]}"
            )

        # Split the report into a tiny metadata JSON + per-table CSVs.
        # The two big arrays (Issue_Summary, Issue_Details) become CSV files that JS reads
        # via DG.DataFrame.fromCsv — much faster than JSON.parse + DG.DataFrame.fromObjects
        # on a 70+ MB JSON blob. Nested arrays in cells (e.g. `variables`, `values`) are
        # JSON-encoded so JS can JSON.parse them per-row when needed.
        appdata_dir = data_path.rstrip("/")
        appdata_json_target = f"{appdata_dir}/validation_results.json"
        appdata_summary_csv_target = f"{appdata_dir}/issue_summary.csv"
        appdata_details_csv_target = f"{appdata_dir}/issue_details.csv"
        try:
            with open(result_file, "r", encoding="utf-8") as f:
                report = json.load(f)
            issue_summary = report.pop("Issue_Summary", []) or []
            issue_details = report.pop("Issue_Details", []) or []
            # The remaining report keeps Conformance_Details, Dataset_Details, Rules_Report
            # — small (< 100 KB typically). JS uses it for headers and lookups.
            slim_json_path = result_file.with_name(result_file.stem + ".slim.json")
            with open(slim_json_path, "w", encoding="utf-8") as f:
                json.dump(report, f, separators=(",", ":"))

            summary_csv_path = result_file.with_name(result_file.stem + ".issue_summary.csv")
            details_csv_path = result_file.with_name(result_file.stem + ".issue_details.csv")
            _write_records_as_csv(issue_summary, summary_csv_path)
            _write_records_as_csv(issue_details, details_csv_path)

            logging.info(
                "Uploading slim metadata (%d bytes), summary CSV (%d bytes, %d rows), details CSV (%d bytes, %d rows)",
                slim_json_path.stat().st_size,
                summary_csv_path.stat().st_size, len(issue_summary),
                details_csv_path.stat().st_size, len(issue_details),
            )
            api.files.upload("System:AppData", appdata_json_target, str(slim_json_path))
            api.files.upload("System:AppData", appdata_summary_csv_target, str(summary_csv_path))
            api.files.upload("System:AppData", appdata_details_csv_target, str(details_csv_path))
        except Exception as e:
            logging.error("Failed to split/upload validation report: %s", e)
            # Fall back to uploading the original JSON so the caller still gets something.
            api.files.upload("System:AppData", appdata_json_target, str(result_file))
            raise
        return appdata_json_target

    finally:
        # Restore .sync_meta.json file if it was moved
        if sync_meta_backup and sync_meta_backup.exists():
            try:
                sync_meta_backup.rename(sync_meta_file)
                logging.info(f"Restored .sync_meta.json to {local_data_path}")
            except OSError as e:
                logging.warning(f"Failed to restore .sync_meta.json: {e}")

        # Clean up the temporary XPT conversion directory, if created
        if converted_dir and converted_dir.exists():
            try:
                shutil.rmtree(converted_dir)
                logging.info(f"Removed temporary XPT directory {converted_dir}")
            except OSError as e:
                logging.warning(f"Failed to remove temporary XPT directory {converted_dir}: {e}")


#name: Check Health
#output: string result
@app.task(name="check_health", base=DatagrokTask)
def check_health() -> str:
    logging.info("check_health invoked")
    proc = _run_core(["--help"], timeout=60)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "Unable to execute CORE binary.")
    return json.dumps({"status": "ok"})


#name: List CDISC CORE Rule Sets
#output: string result
@app.task(name="list_rule_sets", base=DatagrokTask)
def list_rule_sets() -> str:
    """Diagnostic: list the rule sets (standard + version combinations) cached in CORE."""
    logging.info("list_rule_sets invoked")
    proc = _run_core(["list-rule-sets"], timeout=120)
    return json.dumps({
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    })


#name: List CDISC CORE Rules
#input: string standard
#input: string version
#output: string result
@app.task(name="list_rules", bind=True, base=DatagrokTask)
def list_rules(self, standard: str, version: str, **kwargs) -> str:
    """Diagnostic: list the rules CORE will run for a given standard + version."""
    logging.info(f"list_rules invoked: standard={standard}, version={version}")
    if not standard:
        raise RuntimeError("Standard is required.")
    if not version:
        raise RuntimeError("Version is required.")
    proc = _run_core(["list-rules", "-s", standard, "-v", version], timeout=300)
    return json.dumps({
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    })


#name: List CORE Cache Contents
#output: string result
@app.task(name="list_cache_contents", base=DatagrokTask)
def list_cache_contents() -> str:
    """Diagnostic: list files in /app and any obvious cache directories,
    so we can tell whether the CORE github release ships rules data."""
    logging.info("list_cache_contents invoked")
    result: Dict[str, object] = {}
    for candidate in ("/app", "/app/resources", "/app/cache", "/app/_internal"):
        p = Path(candidate)
        if not p.exists():
            result[candidate] = "<does not exist>"
            continue
        try:
            entries = []
            for item in sorted(p.iterdir()):
                kind = "dir" if item.is_dir() else "file"
                size = item.stat().st_size if item.is_file() else None
                entries.append({"name": item.name, "type": kind, "size": size})
            result[candidate] = entries
        except OSError as e:
            result[candidate] = f"<error: {e}>"
    return json.dumps(result, indent=2)

