
rmdir /S /Q build
del .\collision\Collision.pyd


REM CALL "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\amd64\vcvars64.bat"
set PY_VER=3.5
set ARCH=64


set BUILD_CONFIG=Release
set CENV=C:\Users\Eelco\Miniconda2\envs\collision
set SP_DIR=%CENV%\Lib\site-packages

set GENERATOR_NAME=Visual Studio 14 2015

REM pick architecture
if %ARCH%==64 (
	set GENERATOR_NAME=%GENERATOR_NAME% Win64
)

REM tell cmake where Python is
set PYTHON_LIBRARY=%CENV%\libs\python%PY_VER:~0,1%%PY_VER:~2,1%.lib

REM work in build subdir
mkdir build
cd build

cmake ../src -G"%GENERATOR_NAME%" ^
    -Wno-dev ^
    -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% ^
    -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
    -DPYTHON_INCLUDE_DIR:PATH="%CENV%/include" ^
    -DPYTHON_LIBRARY:FILEPATH="%PYTHON_LIBRARY%" ^
	-DEIGEN_INCLUDE_DIR:PATH="%CENV%/Library/include" ^
	-DNUMPY_INCLUDE_DIR:PATH="%SP_DIR%/numpy/core/include" ^
    -DBOOST_ROOT:PATH="%CENV%/Library"

cd..

cmake --build ./build --clean-first --target ALL_BUILD --config %BUILD_CONFIG%

if %errorlevel% neq 0 exit /b %errorlevel%



copy .\build\release\Collision.pyd .\collision

python .\collision\test.py
