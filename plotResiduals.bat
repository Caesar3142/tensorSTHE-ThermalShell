@echo off
setlocal
cd /d "%~dp0"
if exist residuals.html (
  start "" "%~dp0residuals.html"
) else (
  echo residuals.html not found.
  echo From Git Bash / OpenFOAM:  ./plotResiduals --open
)
