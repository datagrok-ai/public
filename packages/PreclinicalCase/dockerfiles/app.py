import json
import logging
import os
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

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


def _sync_directory(remote_path: str, token: str, dest: Path) -> None:
    if not token:
        raise RuntimeError("USER_API_KEY was not provided by Datagrok.")
    if not settings.api_url:
        raise RuntimeError("DATAGROK_API_URL is not configured for the container.")

    logging.info(f"remote_path {remote_path}, token: {token}, dest: {dest}, settings.api_url: {settings.api_url}")
    api = DatagrokClient(api_key=token, base_url=settings.api_url)
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
    local_dir = DATASETS_ROOT / data_path
    logging.info(f"local_dir {local_dir}")
    lock_file = LOCKS_ROOT / f"{data_path.replace('/', '_')}.lock"
    logging.info(f"lock_file {lock_file}")
    with FileLock(str(lock_file), timeout=240):
        _sync_directory(data_path, token, local_dir)
    local_data_path = str(local_dir)

    # Temporarily move .sync_meta.json file if it exists (created by Datagrok file sync)
    # Validation fails if folder contains JSON files that aren't SEND datasets
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
    # Validation fails if folder contains JSON files that aren't SEND datasets
    study_config_file = Path(local_data_path) / "studyConfig.json"
    if study_config_file.exists():
        try:
            study_config_file.unlink()
            logging.info(f"Removed studyConfig.json from {local_data_path}")
        except OSError as e:
            logging.warning(f"Failed to remove studyConfig.json: {e}")

    try:
        # Parse options JSON string if provided
        options_dict: Optional[Dict[str, object]] = None
        if options:
            try:
                options_dict = json.loads(options)
            except json.JSONDecodeError as e:
                raise RuntimeError(f"Invalid JSON in options parameter: {e}")

        logging.info(f"local_data_path {local_data_path}")
        args = _build_validation_args(standard, version, local_data_path, output_format, options_dict)

        self.update_state(meta={"description": "Running validation inside container."})
        proc = _run_core(args)

        if proc.returncode not in (0, 1):
            raise RuntimeError(
                f"core validate failed with exit code {proc.returncode}. stderr: {proc.stderr.strip()}"
            )

        # proc.stdout contains the path to the result JSON file in format: "Output: CORE-Report-2025-12-01T18-51-22"
        logging.info(f"proc.stdout {proc.stdout}")
        stdout_text = proc.stdout.strip()
        if not stdout_text:
            raise RuntimeError("Validation completed but no result file path in stdout")

        # Extract filename after "Output: "
        if not stdout_text.startswith("Output: "):
            raise RuntimeError(f"Unexpected stdout format: {stdout_text}")
        
        result_file_name = stdout_text[len("Output: "):].strip()
        if not result_file_name:
            raise RuntimeError("No filename found after 'Output: ' in stdout")

        # Result file is in the working directory (/app) - use absolute path
        result_file = Path(f"/app/{result_file_name}.json")
        if not result_file.exists():
            raise RuntimeError(f"Result file not found: {result_file}")

        # Read and return the contents of the result file
        with open(result_file, "r", encoding="utf-8") as f:
            result_content = f.read()

        return result_content

    finally:
        # Restore .sync_meta.json file if it was moved
        if sync_meta_backup and sync_meta_backup.exists():
            try:
                sync_meta_backup.rename(sync_meta_file)
                logging.info(f"Restored .sync_meta.json to {local_data_path}")
            except OSError as e:
                logging.warning(f"Failed to restore .sync_meta.json: {e}")


#name: List CDISC CORE Rules
#input: string standard
#input: string version
#output: string result
@app.task(name="list_core_rules", bind=True, base=DatagrokTask)
def list_core_rules(
    self,
    standard: str,
    version: str,
    **kwargs,
) -> str:
    """List CDISC CORE rules for a given standard and version."""
    self.update_state(meta={"description": "Listing CDISC CORE rules."})
    
    if not standard:
        raise RuntimeError("Standard is required.")
    if not version:
        raise RuntimeError("Version is required.")
    
    args: List[str] = ["list-rules", "-s", standard, "-v", version]
    
    self.update_state(meta={"description": "Executing list-rules command."})
    proc = _run_core(args, timeout=300)
    
    if proc.returncode != 0:
        raise RuntimeError(
            f"list-rules failed with exit code {proc.returncode}. stderr: {proc.stderr.strip()}"
        )
    
    # Parse the output - list-rules may return table format or JSON
    output = proc.stdout.strip()
    if not output:
        return json.dumps({"rules": []})
    
    # Try to parse as JSON first
    try:
        parsed = json.loads(output)
        return json.dumps(parsed, indent=2)
    except json.JSONDecodeError:
        # If not JSON, parse table format and convert to structured JSON
        # Table format typically has headers and data rows
        lines = [line.strip() for line in output.split("\n") if line.strip()]
        if not lines:
            return json.dumps({"rules": []})
        
        # Try to detect table structure (header row, separator, data rows)
        rules = []
        headers = None
        data_start_idx = 0
        
        for i, line in enumerate(lines):
            # Look for header row (usually contains column names)
            if i == 0 and ("Rule ID" in line or "rule" in line.lower() or "id" in line.lower()):
                # Try to parse header
                headers = [h.strip() for h in line.split("|") if h.strip()] if "|" in line else line.split()
                data_start_idx = i + 1
                # Skip separator line if present
                if i + 1 < len(lines) and all(c in "-" for c in lines[i + 1].strip()):
                    data_start_idx = i + 2
                break
        
        # Parse data rows
        for line in lines[data_start_idx:]:
            if not line.strip() or all(c in "-" for c in line.strip()):
                continue
            # Parse row data
            if "|" in line:
                values = [v.strip() for v in line.split("|") if v.strip()]
            else:
                values = line.split()
            
            if values:
                # Create rule object from parsed values
                rule_obj: Dict[str, str] = {}
                if headers and len(headers) == len(values):
                    for j, header in enumerate(headers):
                        rule_obj[header] = values[j] if j < len(values) else ""
                else:
                    # Fallback: use positional values
                    rule_obj = {
                        "rule_id": values[0] if len(values) > 0 else "",
                        "description": " ".join(values[1:]) if len(values) > 1 else "",
                    }
                rules.append(rule_obj)
        
        return json.dumps({"rules": rules}, indent=2)


@app.task(name="check_health", base=DatagrokTask)
def check_health() -> str:
    proc = _run_core(["--version"], timeout=60)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "Unable to execute CORE binary.")
    return json.dumps({"status": "ok", "version": proc.stdout.strip()})

