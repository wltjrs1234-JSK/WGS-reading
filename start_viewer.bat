@echo off
cd /d "%~dp0"
if not exist "dist\index.html" (
  echo [오류] dist\index.html 파일이 없습니다. 먼저 npm run build를 실행해주세요.
  pause
  exit /b 1
)
start "" wscript.exe //nologo "%~dp0start_viewer.vbs"
exit /b 0

