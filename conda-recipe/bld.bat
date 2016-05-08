
CALL "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\vcvarsamd64.bat"

rmdir /S /Q build

set BUILD_CONFIG=Release
set PY_VER=3.5

REM pick generator based on python version
if %PY_VER%==2.7 (
    set GENERATOR_NAME=Visual Studio 9 2008
)
if %PY_VER%==3.4 (
    set GENERATOR_NAME=Visual Studio 10 2010
)
if %PY_VER%==3.5 (
    set GENERATOR_NAME=Visual Studio 14 2015
)

REM pick architecture
if %ARCH%==64 (
	set GENERATOR_NAME=%GENERATOR_NAME% Win64
)

REM tell cmake where Python is
set PYTHON_LIBRARY=%PREFIX%\libs\python%PY_VER:~0,1%%PY_VER:~2,1%.lib


REM work in build subdir
mkdir build
cd build

cmake ../src -G"%GENERATOR_NAME%" ^
    -Wno-dev ^
    -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% ^
    -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
    -DPYTHON_INCLUDE_DIR:PATH="%PREFIX%/include" ^
    -DPYTHON_LIBRARY:FILEPATH="%PYTHON_LIBRARY%" ^
	-DBOOST_ROOT="%PREFIX%" ^
	-DEIGEN_INCLUDE_DIR:PATH="%PREFIX%/Library/include" ^
	-DNUMPY_INCLUDE_DIR:PATH="%SP_DIR%/numpy/core/include"

cmake --build . --clean-first --target ALL_BUILD --config %BUILD_CONFIG%
REM cmake --build . --clean-first --target INSTALL --config %BUILD_CONFIG%

if errorlevel 1 exit 1



