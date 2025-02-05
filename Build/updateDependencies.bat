@echo OFF
setlocal enabledelayedexpansion
set XBMP_TOOLS_PATH="%cd%"

rem --------------------------------------------------------------------------------------------------------
rem Set the color of the terminal to blue with yellow text
rem --------------------------------------------------------------------------------------------------------
COLOR 8E
powershell write-host -fore White ------------------------------------------------------------------------------------------------------
powershell write-host -fore Cyan Welcome I am your XBMP_TOOLS dependency updater bot, let me get to work...
powershell write-host -fore White ------------------------------------------------------------------------------------------------------
echo.

:DOWNLOAD_DEPENDENCIES
powershell write-host -fore White ------------------------------------------------------------------------------------------------------
powershell write-host -fore White XBMP_TOOLS - DOWNLOADING DEPENDENCIES
powershell write-host -fore White ------------------------------------------------------------------------------------------------------
echo.

rem ------------------------------------------------------------
rem TINY DDS LOADER
rem ------------------------------------------------------------
rmdir "../dependencies/tinyddsloader" /S /Q
git clone https://github.com/benikabocha/tinyddsloader.git "../dependencies/tinyddsloader"
if %ERRORLEVEL% GEQ 1 goto :PAUSE

rem ------------------------------------------------------------
rem STB IMAGE LOADER
rem ------------------------------------------------------------
rmdir "../dependencies/stb" /S /Q
git clone --recurse-submodules -j8  https://github.com/nothings/stb.git "../dependencies/stb"
if %ERRORLEVEL% GEQ 1 goto :PAUSE

rem ------------------------------------------------------------
rem EXR IMAGE LOADER
rem ------------------------------------------------------------
rmdir "../dependencies/tinyexr" /S /Q
git clone https://github.com/syoyo/tinyexr.git "../dependencies/tinyexr"
if %ERRORLEVEL% GEQ 1 goto :PAUSE


rem ------------------------------------------------------------
rem Finishing
rem ------------------------------------------------------------
:DONE
powershell write-host -fore White ------------------------------------------------------------------------------------------------------
powershell write-host -fore White XBMP_TOOLS - DONE!!
powershell write-host -fore White ------------------------------------------------------------------------------------------------------

:PAUSE
rem if no one give us any parameters then we will pause it at the end, else we are assuming that another batch file called us
if %1.==. pause


