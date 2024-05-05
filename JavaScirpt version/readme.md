# Dynamic-Color-JavaScript

This is a frontend code implementation library based on **d3.js**. We place related dependencies in the same folder, with the main file being **dynamic_color.js**, which can be used as a module like the case in that file.

This implementation is relatively simple, involving some parameter simplification and simpler implementation because JavaScript does not support multithreading like OpenMP. This leads to significantly longer execution times for completely consistent code. Additionally, regarding the spatial distribution, it is not directly implemented like the C++ version because the data format used by the users may vary. We have provided data-related interfaces that allow for straightforward implementation and integration.

The JavaScript version provided is a conceptual framework. We recommend using the C++ version for experiments and trials, as it offers better quality and efficiency.

## Code dependency

The core of this project is implemented based on **d3.js**, specifically using version 7. The relevant code is located in the lib folder and can be used and tested under **Node.js**. In addition, we have also utilized some open-source materials for assistance. 

The instructions for using the dependencies are as follows:

- D3.js: https://d3js.org/ Advanced visualization library for the representation and transformation of colors.:
- static-kdtree: https://github.com/mikolalysenko/static-kdtree We use a k-d tree to quickly find nearest neighbors, which is used in the project to make the name difference more accurate.

The latter needs to be configured in **Node.js**, and we recommend using npm:
```bash
npm install static-kdtree
```

Our code also learned from the work of predecessors, and we are very grateful for their outstanding achievements:
- Colorgorical: https://github.com/connorgr/colorgorical The work provided inspiration and we utilized the part concerning name difference.
- Palettailor: https://github.com/IAMkecheng/palettailor-library It provided us with some inspiration in the use of simulated annealing and name distance.

## Code use
The use of the code can refer to the main function in [dynamic_color.js](https://github.com/Dynamic-Color/Dynamic-Color/blob/main/JavaScirpt%20version/dynamic_color.js), where you create a DynamicColor class and then call the run function.

```javascript
function main() {
    // use case for testing
    let dc = new DynamicColor();
    let basic_colors = [[255, 0, 0], [0, 255, 0], [0, 0, 255]];
    let palette_sizes = [5, 5, 5];
    let modify_pcolors = false;
    /* Function run:
    - basic_colors: [[r, g, b], ...] or [] if need to generate basic colors
    - palette_sizes: [size1, size2...], which means the number of colors in each cluster
    - use_data: whether to use data fitting
     */
    let res = dc.run(basic_colors, palette_sizes, true);
    //let res = dc.run([], palette_sizes, true);
    debugger;
}
```

Regarding data items, you can set the method of data import and the form of loss calculation, for which we provide direct interfaces for implementation. You can also refer to the C++ version, which includes implementations of three different forms: scatter plot, PCP (Parallel Coordinates Plot), and grid visualization.
```javascript
load_data(data) {
    this.data = data;
    // TODO: implement this function for specific data format
}

evaluate_data(palette) {
    return 0.0;
    // TODO: implement this function for specific data format and mode(such as discriminative mode, harmony mode, etc.)
}
```


