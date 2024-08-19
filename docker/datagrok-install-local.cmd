@echo off

REM This script automates the Datagrok local installation and running
REM To see additional actions, run "datagrok-install-local.cmd help"

setlocal enabledelayedexpansion

set GREEN=[0;32m
set RED=[0;31m
set YELLOW=[0;33m
set RESET=[0m
set timeout=30

set compose_config_name=localhost.docker-compose.yaml
set datagrok_public_repo_url=https://raw.githubusercontent.com/datagrok-ai/public/master/docker/%compose_config_name%
set datagrok_local_url=http://localhost:8080/

set script_name=%~f0
for %%I in ("%script_name%") do set "script_dir=%%~dpI"
set compose_config_path="%script_dir%\%compose_config_name%"
echo %compose_config_path%

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
  exit /b
:stop
  call :datagrok_stop
  exit /b
:purge
  call :datagrok_purge
  exit /b
:install
  call :datagrok_install
  exit /b
:reset
  call :datagrok_reset
  exit /b
:start
  call :datagrok_start
  exit /b
:message
  echo %YELLOW%%~1%RESET%
  exit /b

:error
  echo %RED%!!!%~1!!!%RESET%
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
  if not exist %compose_config_path% (
    call :message "Downloading Datagrok config file"
    curl -o %compose_config_path% %datagrok_public_repo_url%
  )
  call :message "Pulling Datagrok images (this can take a while depending on your Internet connection speed)"
  docker compose -f %compose_config_path% --profile all pull
  exit /b

:datagrok_start
  call :check_installation
  if errorlevel 1 (
    call :datagrok_install
  )
  call :message "Starting Datagrok containers"
  docker compose -f %compose_config_path% --project-name datagrok --profile all up -d
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
    docker rmi $(docker images -q datagrok/*)
  )
  exit /b

:help
  echo usage: %script_name% install^|start^|stop^|reset^|purge
  exit /b

:check_docker_daemon
  docker info > nul 2>&1
  if %errorlevel% neq 0 (
    exit /b 1
  )
  exit /b 0
