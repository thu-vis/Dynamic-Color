# Dynamic-Color-C++/Python

This is a pybind11 project, which uses C++ to generate libraries usable in Python. Of course, you can also refer to the invocation method in dynamic_color.cpp to rewrite a main.cpp and compile it. You can modify the cmake file to directly use it in C++.

### Code dependency

This project utilizes cmake to generate pybind11 programs, which can be directly called in Python for C code execution; please clone the code from the corresponding folder of pybind11: https://github.com/pybind/pybind11

Environment requirements: g++ cmake
Open-source library dependencies used in the code: eigen, color-util, nanoflann, nlohmann

Most of these dependencies are already included in the include folder (in source code form). If you wish to install these dependencies in other forms, please refer to:

- eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page#Documentation for vector calculations
- color-util: https://github.com/yuki-koyama/color-util for color systems and color difference calculation // Fixed a small error and raised an issue to the original project
- nanoflann: https://github.com/jlblancoc/nanoflann for fast nearest neighbor search using kdtree
- nlohmann/json: https://github.com/nlohmann/json for JSON reading and storage using C++

### Code start

```
mkdir build
cd build
cmake ..
make
```

Then you will find the packaged library in the build directory, for example, under python version 3.8 in linux as dynamic-color.cpython-38-x86_64-linux-gnu.so

You can use it like a library in Python under version 3.8.

**Very important: Please place the "static" folder in the root directory of the calling environment to generate results quickly.**
You can download it from Google Drive:
https://drive.google.com/drive/folders/1XtNUwa650u4l0qPGOCvZ6DbmdTABvz8M?usp=drive_link.
