import os
import sys

if __name__ == "__main__":
    # dir = input("Please input the directory:")
    sample = input("Please input the string:")
    # 获取本目录下所有的文件名
    old_names = os.listdir()
    # 遍历目录下的文件名
    for old_name in old_names:
        # 跳过本脚本文件
        if old_name != sys.argv[0]:
            name =  old_name.split('_bfb')
            # 用新的文件名替换旧的文件名
            os.rename(old_name,name[0]+sample+name[1])