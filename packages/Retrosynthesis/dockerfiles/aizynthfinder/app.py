import logging
import json
from pathlib import Path
import gc

from aizynthfinder.aizynthfinder import AiZynthFinder
from filelock import FileLock
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings
from datagrok_api import DatagrokClient


settings = Settings(log_level=logging.DEBUG)
app = Celery(settings.celery_name, broker=settings.broker_url)

# TODO: Remove it when grok_spawner is updated
settings.api_url = 'http://host.docker.internal:8082'


CONFIGS_DIR = Path.cwd() / "configs"
CONFIGS_DIR.mkdir(parents=True, exist_ok=True)

DEFAULT_CONFIG = str(CONFIGS_DIR / "default" / "config.yml")

LOCK_DIR = Path("/tmp/sync_locks")
LOCK_DIR.mkdir(parents=True, exist_ok=True)

#name: CalculateRetroSynthesisPaths
#meta.cache: all
#meta.cache.invalidateOn: 0 0 * * 1
#input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
#input: string config
#output: string paths
@app.task(name='run_aizynthfind', bind=True, base=DatagrokTask)
def run_aizynthfind(self, molecule: str, config: str, **kwargs):
    if config == "":
        config_path = DEFAULT_CONFIG
    else:
        self.update_state(meta={"description": "Config synchronization in progress..."})
        config_path = _sync_user_config(config, kwargs.get("USER_API_KEY", None))
        self.update_state(meta={"description": "Calculating retrosynthesis paths..."})

    finder = None
    try:    
        finder = AiZynthFinder(configfile=config_path)
        _silence_aizynth_logs()
        finder.stock.select(finder.stock.items)
        finder.expansion_policy.select_all()
        finder.target_smiles = molecule
        finder.tree_search()
        finder.build_routes()
        finder.routes.compute_scores(*finder.scorers.objects())
        routes_dict = finder.routes.dict_with_extra(include_metadata=True, include_scores=True)
        return json.dumps(routes_dict)
    finally:
        if finder is not None:    
            del finder
            gc.collect()

#name: SyncConfig
#input: string config
@app.task(name='sync_config', base=DatagrokTask)
def sync_config(config: str, **kwargs):
    if config == "":
        return
    _sync_user_config(config, kwargs.get("USER_API_KEY", None))


def _sync_user_config(remote_config_path: str, token: str) -> str:
    if not token:
        raise RuntimeError("Api Token wasn't passed by Datagrok server.")

    if not settings.api_url:
        raise RuntimeError(
            "DATAGROK_API_URL wasn't passed as env variable. Check your version of grok_spawner or manually provide it for the container."
        )

    config_path = CONFIGS_DIR / remote_config_path
    lock_path = LOCK_DIR / (remote_config_path.replace("/", "_") + ".lock")

    with FileLock(str(lock_path), timeout=240):
        api = DatagrokClient(token, settings.api_url)
        api.sync_dir("System:AppData", "Retrosynthesis/configs/" + remote_config_path, str(config_path), recursive=False)
        logging.info(f"Directory {remote_config_path} synchronized successfully with {config_path}")

    config_file = config_path / "config.yml"
    if config_file.exists():
        return str(config_file)
    raise RuntimeError("Configs folder should contain config.yml file")


def _silence_aizynth_logs():
    aiz_logger = logging.getLogger("aizynthfinder")
    aiz_logger.handlers.clear()
    aiz_logger.setLevel(logging.CRITICAL + 1)
    aiz_logger.propagate = False
    aiz_logger.addHandler(logging.NullHandler())
