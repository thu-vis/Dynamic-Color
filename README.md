# Dynamic Color Assignment for Hierarchical Data

https://github.com/Dynamic-Color/Dynamic-Color/assets/165145493/99698358-9e72-4df0-be9c-83a5692b2911

=======================================

Codes for color assignment algorithm described in our paper ["Dynamic Color Assignment for Hierarchical Data"](https://xxxx).

## Code Implementation
We provide two versions of code implementation:
- C++/Python version, written in C++ and supported for Python invocation via pybind11. We recommend this version as it allows for acceleration using OpenMP to quickly obtain the desired colors.
You can refer to the [C++ version](https://github.com/Dynamic-Color/Dynamic-Color/tree/main/C%2B%2B%20version) folder for more details.
- JavaScript version, written in JavaScript and based on d3.js for color generation support, is more suitable for direct use in front-end scenarios.
You can refer to the [JavaScript version](https://github.com/Dynamic-Color/Dynamic-Color/tree/main/JavaScirpt%20version) folder for more details.

## Note
Tested on Linux with Python 3.8. 
- C++/Python version is implemented based on CMake to facilitate cross-platform migration. You can modify the CMakeLists.txt according to the corresponding operating system for compilation. 
- JavaScript version is implemented based on Node.js. Please update the Node.js version if it cannot be used.

## Contact
If you have any problem with our code, feel free to contact
- colordynamic0@gmail.com
or describe your problem in Issues.
