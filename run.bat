@echo off

rem 定义 Visual Studio 2019 相关路径
set "VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\Common7\IDE"
set "VCINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise"
set "FrameworkDir=C:\Windows\Microsoft.NET\Framework"
set "FrameworkVersion=v4.8"
set "FrameworkSDKDir=C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um"

rem 检查 VSINSTALLDIR 是否设置
if "%VSINSTALLDIR%"=="" (
    echo Error: VSINSTALLDIR variable is not set.
    echo SYNTAX: %0
    pause
    exit /b 1
)

rem 如果 VCINSTALLDIR 未设置，使用 VSINSTALLDIR
if "%VCINSTALLDIR%"=="" set "VCINSTALLDIR=%VSINSTALLDIR%"

rem 设置开发环境目录
set "DevEnvDir=%VSINSTALLDIR%"

rem 设置 Visual C++ 工具目录
set "MSVCDir=%VCINSTALLDIR%\VC\Tools\MSVC\14.29.30133"

rem 输出环境设置信息
echo Setting environment for using Microsoft Visual Studio 2019 tools.
echo (If you have another version of Visual Studio or Visual C++ installed and wish
echo to use its tools from the command line, run vcvarsall.bat for that version.)

rem 设置环境变量
set "PATH=%DevEnvDir%;%MSVCDir%\bin\Hostx64\x64;%MSVCDir%\bin;%VCINSTALLDIR%\Common7\Tools;%VCINSTALLDIR%\Common7\Tools\bin;%FrameworkSDKDir%\bin;%FrameworkDir%\%FrameworkVersion%;%PATH%"
set "INCLUDE=%MSVCDir%\include;%MSVCDir%\ATLMFC\INCLUDE;%MSVCDir%\PlatformSDK\include;%FrameworkSDKDir%\include;%INCLUDE%"
set "LIB=%MSVCDir%\ATLMFC\LIB;%MSVCDir%\LIB;%MSVCDir%\PlatformSDK\lib;%FrameworkSDKDir%\lib;%LIB%"

rem 输出 INCLUDE 环境变量，用于调试
echo INCLUDE: %INCLUDE%

rem 切换到项目目录
set "PROJECT_DIR=C:\Users\18858\Documents\_ETART"
cd /d "%PROJECT_DIR%" 2>nul
if %errorlevel% neq 0 (
    echo Error: The specified directory "%PROJECT_DIR%" was not found.
    pause
    exit /b 1
)

rem 调用 nmake 进行编译
set "NMAKE_PATH=%MSVCDir%\bin\Hostx64\x64\nmake.exe"
"%NMAKE_PATH%"
if %errorlevel% neq 0 (
    echo Error: Compilation failed.
    pause
    exit /b 1
)

rem 清理目标文件
del *.obj 2>nul

rem 运行生成的可执行文件
if exist etart.exe (
    etart.exe
) else (
    echo etart.exe not found.
    pause
)