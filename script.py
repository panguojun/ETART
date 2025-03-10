import subprocess
import os
from PIL import Image, ImageTk
import tkinter as tk

# 定义 C++ 代码文件和编译后的可执行文件的名称
cpp_file = "C:/Users/18858/Documents/_ETART/etart.cc"
executable = "generate_bmp.exe" if os.name == "nt" else "generate_bmp"
output_bmp = "C:/Users/18858/Documents/_ETART/output/images0.bmp"

# 编译 C++ 代码
try:
    subprocess.run(["g++", cpp_file, "-o", executable], check=True)
    print("C++ 代码编译成功！")
except subprocess.CalledProcessError as e:
    print(f"编译失败: {e}")
    exit(1)

# 运行编译后的可执行文件
try:
    subprocess.run([executable], check=True)
    print("C++ 程序运行成功，已生成 BMP 图像！")
except subprocess.CalledProcessError as e:
    print(f"运行失败: {e}")
    exit(1)

# 检查 BMP 文件是否存在
if not os.path.exists(output_bmp):
    print("未找到生成的 BMP 图像文件。")
    exit(1)

# 使用 Pillow 打开 BMP 图像
try:
    image = Image.open(output_bmp)
except Exception as e:
    print(f"无法打开 BMP 图像: {e}")
    exit(1)

# 创建 Tkinter 窗口
root = tk.Tk()
root.title("显示生成的 BMP 图像")

# 将 Pillow 图像转换为 Tkinter 可用的图像对象
tk_image = ImageTk.PhotoImage(image)

# 创建一个标签并显示图像
label = tk.Label(root, image=tk_image)
label.pack()

# 运行 Tkinter 主循环
root.mainloop()