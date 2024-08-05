# Dynamic Color Assignment for Hierarchical Data

https://github.com/Dynamic-Color/Dynamic-Color/assets/165145493/99698358-9e72-4df0-be9c-83a5692b2911

=======================================

Codes for dynamic color assignment described in our paper ["Dynamic Color Assignment for Hierarchical Data"](https://arxiv.org/abs/2407.14742).

## Code Implementation
We provide two versions of code implementation:

- C++/Python version, written in C++ and supported for Python invocation via pybind11. We recommend this version as it allows for acceleration using OpenMP to quickly obtain the desired colors.
You can refer to the [C++ version](https://github.com/Dynamic-Color/Dynamic-Color/tree/main/C%2B%2B%20version) folder for more details.
- JavaScript version, written in JavaScript and based on d3.js for color generation support, is more suitable for direct use in front-end scenarios.
You can refer to the [JavaScript version](https://github.com/Dynamic-Color/Dynamic-Color/tree/main/JavaScirpt%20version) folder for more details.

Based on our code implementation, high-quality color assignment results for hierarchical data can be generated, which simultaneously considers color discriminability, harmony, and the spatial distribution of data.
![图片](https://github.com/Dynamic-Color/Dynamic-Color/assets/165145493/8be030fa-6cc3-43c4-942e-7d04e369c210)

## Tips for usage
It must be acknowledged that **our method involves a degree of randomness and is influenced by the parameters used** (due to the implementation of multi-threading with OpenMP and the simulated annealing algorithm). The quantative experiments presented in the paper are the average results of multiple trials. The examples provided are relatively representative cases that better demonstrate the characteristics of our method: **ours-D increases the color difference between adjacent elements, while ours-S makes similar categories have similar colors**. If you need to use our work in your project, you may need to be aware of the following parameter adjustment recommendations.
- Firstly, the **color range** needs to be adapted to the environment in which it will be used, especially the **lightness**. You can modify this at the beginning of [scope.hpp](https://github.com/Dynamic-Color/Dynamic-Color/blob/main/C%2B%2B%20version/include/scope.hpp). Our recommendations regarding this are:
  - If you intend to use it in a PDF document or on a darker screen, we recommend a higher lightness/chroma range, such as [50, 85],
  - If you intend to use it for data visualization on a brighter screen, we recommend a looser range, such as [40, 85].
- Secondly, the **variable tradeoff** controls the relative importance of between-class and inter-class differences. The default value is 1. You can modify this parameter to suit your specific environment:
  - If you want the colors under different parents to be as distinct as possible, you can increase the tradeoff value to 2.
  - If you want the subclasses to be as distinct as possible while allowing for similar colors under different parents, you can decrease the tradeoff value, for example, to 0.5.
- Finally, although our method exhibits some randomness, conducting **multiple experiments** can help you find more suitable results.If you still encounter difficulties with adapting the parameters, feel free to **contact us**.

## Note
Tested on Linux with Python 3.8. 
- C++/Python version is implemented based on CMake to facilitate cross-platform migration. You can modify the CMakeLists.txt according to the corresponding operating system for compilation. 
- JavaScript version is implemented based on Node.js. Please update the Node.js version if it cannot be used.

## Contact
If you have any problem with our code, feel free to contact
- jiashu0717c@gmail.com

or describe your problem in Issues.
