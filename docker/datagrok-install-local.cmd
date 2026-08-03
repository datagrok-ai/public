@echo off

REM This script automates the Datagrok local installation and running
REM To see additional actions, run "datagrok-install-local.cmd help"

setlocal enabledelayedexpansion

set GREEN=[0;32m
set RED=[0;31m
set YELLOW=[0;33m
set RESET=[0m
set timeout=30

REM Any tag of datagrok/datagrok: `latest` for the newest stable release, an exact `X.Y.Z`,
REM `bleeding-edge` to track master, or an `X.Y.Z-rc` candidate.
set datagrok_default_version=latest

REM A compose file and the images it starts must come from the same release: the compose file
REM carries the GROK_PARAMETERS that the datlas build inside datagrok/datagrok parses, and the
REM two evolve together. So the compose file is not taken from whatever master happens to hold
REM -- it is taken from the branch the requested image was actually built from, which the image
REM records in its `branch` label. That works for every tag, including moving ones like
REM `latest`, and it means the pairing is derived from the image rather than guessed from the
REM tag's name. resolve_release fills in resolved_ref and resolved_release.
if "%DATAGROK_VERSION%"=="" set DATAGROK_VERSION=%datagrok_default_version%

set compose_config_name=localhost.docker-compose.yaml
set datagrok_image=datagrok/datagrok:%DATAGROK_VERSION%
set datagrok_local_url=http://localhost:8080/
set "resolved_ref="
set "resolved_release="
set image_stale=false

set script_name=%~f0
for %%I in ("%script_name%") do set "script_dir=%%~dpI"
set compose_config_path="%script_dir%%compose_config_name%"

call :check_docker_daemon
if %errorlevel% neq 0 (
    call :message "Docker daemon is not running, please launch Docker Desktop application."
    pause
    exit /b
)

if "%~1"=="" goto :run_application

echo Waiting ...

echo.& call :message "This script automates the Datagrok local installation and running"
echo.& call :message "To see additional actions, run 'datagrok-install-local.cmd help'"

set action=%~1
shift

goto :%action%

:run_application
  call :datagrok_start
  exit /b %errorlevel%
:stop
  call :datagrok_stop
  exit /b
:purge
  call :datagrok_purge
  exit /b
:install
  call :datagrok_install
  exit /b %errorlevel%
:reset
  call :datagrok_reset
  exit /b
:start
  call :datagrok_start
  exit /b %errorlevel%
:update
  call :datagrok_update
  exit /b %errorlevel%
:version
  call :resolve_release
  echo %DATAGROK_VERSION% -^> !resolved_release! (!resolved_ref!)
  exit /b
:message
  echo %YELLOW%%~1%RESET%
  exit /b

:error
  echo %RED%ERROR: %~1%RESET%
  exit /b

:user_query_yn
  echo.%YELLOW%
  set /p "answer=%~1 (y/N)"
  echo.%RESET%
  if /i "!answer!"=="y" (
    exit /b 0
  ) else (
    exit /b 1
  )

:count_down
  echo Waiting:
  setlocal enabledelayedexpansion
  for /l %%I in (%1,-1,1) do (
    echo %%I
    timeout 1 > nul
  )
  echo 0
  endlocal
  exit /b

:check_docker
  where docker > nul 2>&1
  if errorlevel 1 (
    call :error "Docker engine is not installed"
    call :message "Please install Docker and Docker Compose plugin manually: https://docs.docker.com/engine/install/"
    exit /b 255
  )
  exit /b

REM What the registry says, into registry_branch. For a moving tag like `latest` this is the
REM only authority -- a copy on disk may point at an older release. imagetools reads the image
REM config straight from the registry, so this costs no download.
:registry_branch
  set "registry_branch="
  for /f "delims=" %%I in ('docker buildx imagetools inspect %datagrok_image% --format "{{ index .Image.Config.Labels \"branch\" }}" 2^>nul') do set "registry_branch=%%I"
  if /i "!registry_branch!"=="<no value>" set "registry_branch="
  exit /b 0

REM What the copy on disk says, into local_branch. Empty if it is absent or predates the label.
:local_branch
  set "local_branch="
  for /f "delims=" %%I in ('docker image inspect %datagrok_image% --format "{{.Config.Labels.branch}}" 2^>nul') do set "local_branch=%%I"
  if /i "!local_branch!"=="<no value>" set "local_branch="
  exit /b 0

REM Work out which branch built the image this tag names. Registry first; a local image is the
REM fallback for when the registry cannot be reached, so an existing install still starts offline.
:resolve_release
  if not "!resolved_ref!"=="" exit /b 0
  call :registry_branch
  call :local_branch
  set "image_branch=!registry_branch!"
  REM A local copy older than the tag now points at has to be replaced: running it against the
  REM newer release's compose file is exactly the pairing this all exists to prevent.
  if not "!registry_branch!"=="" if not "!local_branch!"=="" (
    if /i not "!registry_branch!"=="!local_branch!" set image_stale=true
  )
  if "!image_branch!"=="" set "image_branch=!local_branch!"
  if "!image_branch!"=="" (
    call :message "Resolving %DATAGROK_VERSION%"
    docker pull %datagrok_image% > nul
    if errorlevel 1 (
      call :error "Could not pull %datagrok_image%"
      call :message "Check that '%DATAGROK_VERSION%' is a published tag of datagrok/datagrok."
      exit /b 255
    )
    call :local_branch
    set "image_branch=!local_branch!"
  )
  if "!image_branch!"=="" (
    REM Images built before the label existed: read the tag instead.
    if /i "%DATAGROK_VERSION%"=="bleeding-edge" (
      set "resolved_ref=master"
      set "resolved_release=bleeding-edge"
    ) else (
      for /f "tokens=1 delims=-" %%V in ("%DATAGROK_VERSION%") do (
        set "resolved_ref=release/%%V"
        set "resolved_release=%%V"
      )
    )
    call :message "Note: %datagrok_image% has no branch label; assuming !resolved_release!."
  ) else (
    echo !image_branch! | findstr /b /c:"release/" > nul
    if errorlevel 1 (
      REM A build from master or a feature branch: those track master's compose file.
      set "resolved_ref=master"
      set "resolved_release=bleeding-edge"
    ) else (
      set "resolved_ref=!image_branch!"
      set "resolved_release=!image_branch:release/=!"
    )
  )
  set datagrok_public_repo_url=https://raw.githubusercontent.com/datagrok-ai/public/!resolved_ref!/docker/%compose_config_name%
  if /i not "%DATAGROK_VERSION%"=="!resolved_release!" (
    call :message "%DATAGROK_VERSION% is !resolved_release! (built from !resolved_ref!)"
  )
  exit /b 0

REM Reads the release a compose file on disk was cut from into compose_release.
:compose_release
  set "compose_release="
  if not exist %compose_config_path% exit /b 0
  for /f "tokens=1,* delims=:" %%A in ('findstr /b /c:"x-datagrok-release:" %compose_config_path%') do (
    for /f "tokens=* delims= " %%C in ("%%B") do set "compose_release=%%C"
  )
  exit /b 0

REM Refuse to start a compose file from one release against images from another. The
REM comparison is against the release the image resolved to, not the tag the user typed, so a
REM moving tag like `latest` is checked against what it currently points at. Files cut before
REM the marker existed have nothing to compare; those are only warned about.
:verify_compose_release
  call :compose_release
  if "!compose_release!"=="" (
    call :message "Note: %compose_config_name% predates release binding; assuming !resolved_release!."
    exit /b 0
  )
  if /i not "!compose_release!"=="!resolved_release!" (
    call :error "Version mismatch"
    call :message "%compose_config_name% is built for '!compose_release!', but %DATAGROK_VERSION% is '!resolved_release!'."
    call :message "Running them together makes the Datagrok server fail to start."
    call :message "Delete %compose_config_path% and re-run this script to download the matching file."
    exit /b 255
  )
  exit /b 0

:download_compose
  call :message "Downloading Datagrok config file for !resolved_release!"
  curl -fsSL -o %compose_config_path% !datagrok_public_repo_url!
  if errorlevel 1 (
    call :error "Could not download the config file"
    call :message "Tried: !datagrok_public_repo_url!"
    call :message "%datagrok_image% was built from !resolved_ref!, but that branch has no compose file."
    exit /b 255
  )
  call :verify_compose_release
  exit /b %errorlevel%

:check_installation
  if not exist %compose_config_path% (
    exit /b 1
  ) else (
    for /f "delims=" %%I in ('docker images -q "datagrok/datagrok"') do (
      if not "%%I"=="" (
        exit /b 0
      )
    )
    exit /b 1
  )

:datagrok_install
  call :resolve_release
  if errorlevel 1 exit /b 255
  REM A file left over from another release is replaced, not reused.
  call :compose_release
  if not exist %compose_config_path% (
    call :download_compose
    if errorlevel 1 exit /b 255
  ) else if /i not "!compose_release!"=="!resolved_release!" (
    if not "!compose_release!"=="" (
      call :download_compose
      if errorlevel 1 exit /b 255
    )
  )
  call :verify_compose_release
  if errorlevel 1 exit /b 255
  call :message "Pulling Datagrok !resolved_release! images (this can take a while depending on your Internet connection speed)"
  docker compose -f %compose_config_path% --profile all pull
  exit /b

:datagrok_start
  call :check_installation
  if errorlevel 1 (
    call :datagrok_install
    if errorlevel 1 exit /b 255
  )
  call :resolve_release
  if errorlevel 1 exit /b 255
  call :verify_compose_release
  if errorlevel 1 exit /b 255
  REM The tag moved since this image was pulled; take the new one so the running server
  REM matches the compose file about to configure it.
  if "!image_stale!"=="true" (
    call :message "%DATAGROK_VERSION% has moved to !resolved_release!; pulling"
    docker compose -f %compose_config_path% --profile all pull
  )
  call :message "Starting Datagrok containers"
  docker compose -f %compose_config_path% --project-name datagrok --profile all up -d
  call :report_started
  exit /b

:datagrok_update
  call :resolve_release
  if errorlevel 1 exit /b 255
  call :download_compose
  if errorlevel 1 exit /b 255
  call :message "Updating Datagrok to !resolved_release!"
  docker compose -f %compose_config_path% --project-name datagrok --profile all pull
  docker compose -f %compose_config_path% --project-name datagrok --profile all up -d --force-recreate
  call :report_started
  call :message "Removing old images"
  docker image prune -f
  exit /b

:report_started
  call :message "Waiting while the Datagrok server is starting"
  echo When the browser opens, use the following credentials to log in:
  echo ------------------------------
  echo %GREEN%Login:    admin
  echo Password: admin
  echo %RESET%------------------------------
  echo If you see the message 'Datagrok server is unavailable', just wait for a while and reload the web page
  call :count_down %timeout%
  call :message "Running browser"
  start "" "%datagrok_local_url%"
  call :message "If the browser doesn't open, use the following link: %datagrok_local_url%"
  call :message "To extend Datagrok functionality, install extension packages on the 'Manage -> Packages' page"
  exit /b

:datagrok_stop
  call :check_installation
  if errorlevel 1 (
    call :message "The Datagrok installation was not found. Nothing to stop."
    exit /b 255
  )
  call :message "Stopping Datagrok containers"
  docker compose -f %compose_config_path% --project-name datagrok --profile all stop
  FOR /F "tokens=*" %%i IN ('docker ps --format "{{.Names}}" ^| find "datagrok"') DO docker rm -f %%i
  exit /b

:datagrok_reset
  call :check_installation
  if errorlevel 1 (
    call :message "The Datagrok installation was not found. Nothing to reset."
    exit /b 255
  )
  call :user_query_yn "This action will stop Datagrok, remove all user settings, data, and installed packages. Are you sure?"
  if errorlevel 1 (
    docker compose -f %compose_config_path% --project-name datagrok --profile all stop
    docker compose -f %compose_config_path% --project-name datagrok --profile all down --volumes
  )
  exit /b

:datagrok_purge
  call :check_installation
  if errorlevel 1 (
    call :message "The Datagrok installation was not found. Nothing to remove."
    exit /b 255
  )
  call :user_query_yn "This action will stop Datagrok and COMPLETELY remove the Datagrok installation. Are you sure?"
  if errorlevel 1 (
    docker compose -f %compose_config_path% --project-name datagrok --profile all down --volumes
    for /f "delims=" %%I in ('docker images -q datagrok/*') do docker rmi %%I
  )
  exit /b

:help
  echo usage: %script_name% install^|start^|stop^|update^|reset^|purge^|version
  echo.
  echo   DATAGROK_VERSION is any tag of datagrok/datagrok (default %datagrok_default_version%):
  echo     latest         newest stable release
  echo     ^<X.Y.Z^>        that exact release
  echo     ^<X.Y.Z^>-rc     that release candidate
  echo     bleeding-edge  latest build from master
  echo.
  echo   The compose file is taken from the branch the chosen image was built from,
  echo   so the configuration always matches the images.
  exit /b

:check_docker_daemon
  docker info > nul 2>&1
  if %errorlevel% neq 0 (
    exit /b 1
  )
  exit /b 0
