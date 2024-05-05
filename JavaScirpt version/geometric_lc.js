// geometric lc implementation
// geometric_lc: geometric lightness-chroma harmony
var d3 = require("./lib/d3.v7.min.js");

const kc = 1;
const kl = 2;

var geo_lc_linear = function (colors) {
    let hcl_colors = colors.map(d => d3.hcl(d));
    let _c = 0, _l = 0, _w = 0;
    hcl_colors.forEach((color) => {
        let sc = 1 + 0.045 * color.c;
        let sl = 1 + 0.015 * Math.pow(color.l - 50, 2) / Math.sqrt(20 + Math.pow(color.l - 50, 2));
        let w = Math.pow(kc * sc * kl * sl, -2);
        color.w = w;
        _c += w * color.c;
        _l += w * color.l;
        _w += w;
    });

    _c = _c / _w;
    _l = _l / _w;
    let tp1 = 0, tp2 = 0;
    hcl_colors.forEach((color) => {
        tp1 += color.w * (_l - color.l) * (_c - color.c);
        tp2 += color.w * (Math.pow(_l - color.l, 2) - Math.pow(_c - color.c, 2));
    });
    let _phi = 0.5 * Math.atan2(-2 * tp1, tp2);
    let cos_phi = Math.cos(_phi);
    let sin_phi = Math.sin(_phi);
    let _r = _c * cos_phi + _l * sin_phi;

    let _MD = 0;
    hcl_colors.forEach((color) => {
        color.MD = Math.abs(color.c * cos_phi + color.l * sin_phi - _r);
        if (color.MD > 15.0) _MD += (color.MD - 15.0);
    });
    return _MD;
}

module.exports = geo_lc_linear;
// export default geo_lc_linear;   