# 定义 Visual Studio 2019 相关路径
MSVCDIR = C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.29.30133
WINDOWSKITSDIR = C:\Program Files (x86)\Windows Kits\10

# 定义包含目录，添加 ucrt 目录
INCLUDE_DIRS = /I"$(MSVCDIR)\include" /I"$(MSVCDIR)\ATLMFC\INCLUDE" /I"$(WINDOWSKITSDIR)\Include\10.0.19041.0\um" /I"$(WINDOWSKITSDIR)\Include\10.0.19041.0\shared" /I"$(WINDOWSKITSDIR)\Include\10.0.19041.0\ucrt"

# 定义目标
all: etart.exe

# 生成可执行文件的规则
etart.exe: etart.obj
    @echo "Linking etart.exe..."
    link /OUT:etart.exe etart.obj /LIBPATH:"$(MSVCDIR)\lib\x64" /LIBPATH:"$(WINDOWSKITSDIR)\Lib\10.0.19041.0\um\x64" /LIBPATH:"$(WINDOWSKITSDIR)\Lib\10.0.19041.0\ucrt\x64"

etart.obj: etart.cpp
    @echo "Compiling etart.cpp..."
    cl /Fe:etart.obj $(INCLUDE_DIRS) etart.cpp

# 清理规则，用于删除生成的目标文件和可执行文件
clean:
    @echo "Cleaning up..."
    -del etart.exe etart.obj
    @echo "Cleanup completed."