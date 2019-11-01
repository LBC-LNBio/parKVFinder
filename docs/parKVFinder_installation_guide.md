# KVFinder Installation Guide

## Linux and MacOS

On the terminal, go to the directory containing **KVFinder_linux_distribution.tar.gz**. Then, follow these steps:

1) Untar *KVFinder_linux_distribution.tar.gz*.
```bash
$ tar -zxvf KVFinder_linux_distribution.tar.gz
```
2) Change working directory to KVFinder_linux folder (Note: folder will show after decompression)*.
```bash    
$ cd KVFinder_linux
```
3) Export *KVFinder_PATH* variable and add to *~/.bashrc* file.
```bash
$ export KVFinder_PATH=`pwd`
$ echo "export KVFinder_PATH=`pwd` >> ~/.bashrc
```
4) Recompile KVFinder.
```bash
$ make clean
$ make
```

### PyMOL v1.8 Plug-in installation

A graphical interface is available to use KVFinder alongside PyMOL v1.8.
To install the plug-in, follow these steps:

1) Open PyMOL v1.8
2) Go to **Plugin** menu &rarr; **Plugin Manager**.
3) **Plugin Manager** window will open, go to **Install New Plugin** tab.
4) Under **Install from local file** Group, click in **Choose file...** .
5) **Install Plugin** window will open, go to **KVFinder_linux** directory and select **kvfinder_pymolplugin_linux.py**.
6) **Select plugin directory** window will open, select **/home/\<user\>/.pymol/startup** and click **OK**.
7) **Confirm** window will open, click **OK**.
8) **Sucess** window will appear confirming that the plugin has been installed.
9) Restart PyMOL v1.8.
10) KVFinder plugin is ready to use.
