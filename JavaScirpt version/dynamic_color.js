const d3 = require("./lib/d3.v7.min");
const d3color = require('./lib/d3.color.min');
const kdtree = require('static-kdtree');
const ColorScope = require('./scope');
const fs = require("fs");
const utils = require('./lib/utils');
const matsuda_templates = require('./hue_templates');
const geo_lc_linear = require('./geometric_lc');

// ------------ prepare for color distance--------------------------------
// This section refers to the previous use of namedifference by paletailor and others. 
// FROM: https://github.com/IAMkecheng/palettailor-library Thanks very much.
let c3 = {};
function c3_init(json) {
    var i, C, W, T, A, ccount, tcount;

    // parse colors
    c3.color = [];
    for (i = 0; i < json.color.length; i += 3) {
        c3.color[i / 3] = d3.lab(
            json.color[i],
            json.color[i + 1],
            json.color[i + 2]
        );
    }
    C = c3.color.length;

    // parse terms
    c3.terms = json.terms;
    W = c3.terms.length;

    // parse count table
    c3.T = T = [];
    for (var i = 0; i < json.T.length; i += 2) {
        T[json.T[i]] = json.T[i + 1];
    }

    // construct counts
    c3.color.count = ccount = [];
    for (i = 0; i < C; ++i) ccount[i] = 0;
    c3.terms.count = tcount = [];
    for (i = 0; i < W; ++i) tcount[i] = 0;
    d3.range(T.length).forEach(function (idx) {
        var c = Math.floor(idx / W),
            w = Math.floor(idx % W),
            v = T[idx] || 0;
        ccount[c] += v;
        tcount[w] += v;
    });

    // parse word association matrix
    c3.A = A = json.A;
}
function c3_api() {
    var C = c3.color.length,
        W = c3.terms.length,
        T = c3.T,
        A = c3.A,
        ccount = c3.color.count,
        tcount = c3.terms.count;

    c3.color.cosine = function (a, b) {
        var sa = 0,
            sb = 0,
            sc = 0,
            ta,
            tb;
        for (var w = 0; w < W; ++w) {
            ta = T[a * W + w] || 0;
            tb = T[b * W + w] || 0;
            sa += ta * ta;
            sb += tb * tb;
            sc += ta * tb;
        }
        return sc / Math.sqrt(sa * sb);
    };

    c3.color.vector = function (a) {
        var v = [];
        for (var w = 0; w < W; ++w) {
            v.push(T[a * W + w] || 0);
        }
        return v;
    };
}
c3 = { version: '1.0.0' };
// use kdtree to find color name index
var color_name_map = {};
var kd_points = [];
var kd;
c3.load = function (uri) {
    const filePath = uri;
    fs.readFile(filePath, "utf8", (err, data) => {
        if (err) {
            console.error("Error reading file:", err);
            return;
        }
        c3_init(JSON.parse(data));
        c3_api();
        for (var c = 0; c < c3.color.length; ++c) {
            var x = c3.color[c];
            color_name_map[[x.L, x.a, x.b].join(",")] = c;
        }

        for (var c = 0; c < c3.color.length; ++c) {
            var x = c3.color[c];
            kd_points.push([x.l, x.a, x.b]);
            color_name_map[[x.l, x.a, x.b].join(',')] = c;
        }

        kd = kdtree(kd_points);
        main();
    });
};
c3.load('./lib/c3_data.json');

function getColorNameIndex(c) {
    var x = d3mlab(c),
        index = kd.nn([x.L, x.a, x.b]),
        s = kd_points[index].join(',');
    return color_name_map[s];
}
function getNameDifference(x1, x2) {
    let c1 = getColorNameIndex(x1),
        c2 = getColorNameIndex(x2);
    return 1 - c3.color.cosine(c1, c2);
}

function cdist(x1, x2) {
    return d3color.d3_ciede2000(d3mlab(x1), d3mlab(x2));
}

const max = utils.max;
const min = utils.min;
const fmod = utils.fmod;
const randOne = utils.random;
const randInt = utils.randInt;

// ------------ prepare for color constraints--------------------------------
// modify color fit constraints
function hdist(h1, h2) {
    let d = Math.abs(h1 - h2);
    return Math.min(d, 360 - d);

}
function modifyColorInColorRange(color, base_color, range, hue_min = 0, hue_max = 0) {
    // fit color in range of base_color
    if (range < 1e-5) {
        return d3.rgb(base_color);
    }
    let dis = cdist(color, base_color);
    let ratio = dis > 0 ? range / dis : 1.0;
    let hcl = d3.hcl(color);
    let hcl_base = d3.hcl(base_color);
    let ct = 0;
    while (dis > range && ct <= 20) {
        ct += 1;
        let h = hcl.h + (hcl_base.h - hcl.h) * ratio;
        let c = hcl.c + (hcl_base.c - hcl.c) * ratio;
        let l = hcl.l + (hcl_base.l - hcl.l) * ratio;
        color = d3.hcl(h, c, l);
        dis = cdist(color, base_color);
    }
    if (Math.abs(hue_max - hue_min) > 1e-5) {
        if ((hue_max > hue_min && (hcl.h < hue_min || hcl.h > hue_max)) ||
            (hue_max < hue_min && hcl.h > hue_max && hcl.h < hue_min)) {
            if (hdist(hcl.h, hue_min) < hdist(hcl.h, hue_max)) {
                color = d3.hcl(hue_min, hcl.c, hcl.l);
            }
            else {
                color = d3.hcl(hue_max, hcl.c, hcl.l);
            }
        }
    }
    return color;
}

function modifyColorInConstraint(color, base_id, base_colors) {
    let base_color = base_colors[base_id];
    let dis = cdist(color, base_color);
    let min_dist = dis;
    let min_id = base_id;
    base_colors.forEach((color, i) => {
        let dis = cdist(color, base_color);
        if (dis < min_dist) {
            min_dist = dis;
            min_id = i;
        }
    });
    if (min_id != base_id) {
        let ratio = (dis - min_dist) / dis;
        let hcl = d3.hcl(color);
        let hcl_base = d3.hcl(base_colors[min_id]);
        let h = hcl.h + (hcl_base.h - hcl.h) * ratio;
        let c = hcl.c + (hcl_base.c - hcl.c) * ratio;
        let l = hcl.l + (hcl_base.l - hcl.l) * ratio;
        color = d3.hcl(h, c, l);
    }
    return color;
}

// ------------ prepare for color harmony--------------------------------
const templates = new matsuda_templates();
function template_score(colors) {
    let hues = colors.map(d => d3.hsl(d).h);
    return templates.check_templates(hues)[0] * Math.PI / 180;
}

function set_min_hue(palette, hue_delta) {
    // set min hue delta for parent colors
    let hues = palette.map(d => d3.hcl(d).h);
    let arg_sort = utils.argsort(hues);
    let resort = [];
    arg_sort.forEach((d, i) => {
        resort[d] = i;
    });
    hues.sort((a, b) => a - b);
    let new_hues = [];
    new_hues.push(hues[0]);
    let delta0 = hues[0] + 360 - hues[hues.length - 1];
    for (let i = 1; i < hues.length; i++) {
        let delta = max(hue_delta, hues[i] - hues[i - 1]);
        let new_hue = new_hues[i - 1] + delta;
        new_hues.push(new_hue);
    }
    let delta = (hue_delta > delta0) ? delta0 : hue_delta;
    let new_hue = hues[hues.length - 1] + delta;
    let r = 1.0;
    if (new_hue - hues[0] > 360) {
        r = 360.0 / (new_hue - hues[0]);
    }
    if (r < 1.0) {
        new_hues = new_hues.map(d => d * r);
    }
    // rotate
    let sum_delta = 0;
    new_hues.forEach((d, i) => {
        sum_delta += d - hues[i];
    });
    sum_delta /= hues.length;
    new_hues = new_hues.map(d => {
        let h = d - sum_delta;
        return fmod(h, 360.0);
    });

    let new_palette = [];
    palette.forEach((d, i) => {
        let hcl = d3.hcl(d);
        new_palette.push(d3.hcl(new_hues[resort[i]], hcl.c, hcl.l));
    });
    return new_palette;
}

// ------------ dynamic color generator--------------------------------
class DynamicColor {
    constructor() {
        this.gen_color = false;

        this.scope = new ColorScope();
        this.clusters = [];
        this.idx2id = {};
        this.iter_num = 32;
        // kernel args
        this.global_color_dist = 10;
        this.tradeoff = 1.0;
        this.harg = 1.0;

        // attribute estimate args
        this.estimate_delta = 10.0;

        this.hue_max_range = 180.0;
        this.hue_min_range = 5.0;
        this.init_hue_min_delta = 20.0;

        this.color_max_range = 100.0;
        this.color_min_range = 10.0
        this.init_color_min_delta = 20.0;
        this.shift_p = 1.0;
        this.harmony_args = [1.0, 1.0];
        this.bharg = 0.8; // Most works consider hue harmony to be more important, so we recommend a normalized balance parameter of 0.8. Of course, 1 is also possible

        // add data fitting
        this.data = {};
    }

    /*Functional api:
    - basic_colors: [[r, g, b], ...] or [] if need to generate basic colors
    - palette_sizes: [size1, size2...], which means the number of colors in each cluster
    - use_data: whether to use data fitting
     */
    run(basic_colors, palette_sizes, use_data = false) {
        let modified_basic_colors = this.findBasicColor(basic_colors, palette_sizes, use_data);
        modified_basic_colors = modified_basic_colors.map(d => {
            let color = d3.rgb(d);
            return [color.r, color.g, color.b];
        });
        return this.findColor(modified_basic_colors, palette_sizes, false, use_data);
    }


    // findbasicColor: find parent colors
    findBasicColor(basic_colors, palette_sizes, use_data = false) {
        if (basic_colors.length == 1) {
            this.gen_color = false;
            return [d3.rgb(basic_colors[0][0], basic_colors[0][1], basic_colors[0][2])];
        }
        this.scope.setMode(true);
        let res = [];
        this.gen_color = false;
        if (basic_colors.length == 0 && palette_sizes.length > 0) {
            if (palette_sizes.length == 1) {
                res.push(d3.rgb(254, 228, 179));
            }
            else {
                let init_rgb = d3.rgb(d3.hcl(this.scope.randomHue(), 65, 65));
                let tp_palette = [[init_rgb.r, init_rgb.g, init_rgb.b]];
                let tp_size = [palette_sizes.length];
                res = this.findColor(tp_palette, tp_size, true, use_data);
            }
            this.gen_color = true;
        }
        else {
            let tp_size = palette_sizes.map(d => 1);
            res = set_min_hue(basic_colors.map(d => d3.rgb(d[0], d[1], d[2])), this.init_hue_min_delta);
            res = this.findColor(res.map(d => {
                let rgb = d3.rgb(d);
                return [rgb.r, rgb.g, rgb.b];
            }), tp_size, true, use_data);
            this.gen_color = false;
        }
        this.scope.setMode(false);
        return res;
    }

    // findcolor: unchange input basic colors
    findColor(basic_colors, palette_sizes, modify_pcolors = false, use_data = false) {
        // 1. init color and color range 
        basic_colors = basic_colors.map(color => this.scope.modifyColorInScope(d3.rgb(color[0], color[1], color[2])));
        this.clusters = [];
        let start = 0, end = 0;
        let all_hues = [], sqrt_nums = [];
        basic_colors.forEach((color, index) => {
            start = end;
            end = start + palette_sizes[index];
            this.clusters.push({
                color: color,
                num: palette_sizes[index],
                start: start,
                end: end,
                index: index
            });
            for (let i = start; i < end; i++) {
                this.idx2id[i] = index;
            }
            all_hues.push(color.h);
            sqrt_nums.push(Math.sqrt(palette_sizes[index]));
        });
        if (all_hues.length > 1) {
            // get argsort of all_hues
            let arg_hues = utils.argsort(all_hues);
            let length = all_hues.length;
            let fi = length - 1;
            for (let i = 0; i < length; i++) {
                // by neighbor delta hue
                let min_range;
                if (i == 0) {
                    let gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[fi]]) * this.tradeoff;
                    min_range = (360.0 - all_hues[arg_hues[fi]] + all_hues[arg_hues[i]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[fi]]));
                }
                else {
                    let gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[i - 1]]) * this.tradeoff;
                    min_range = (all_hues[arg_hues[i]] - all_hues[arg_hues[i - 1]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[i - 1]]));
                }

                let max_range;
                if (i == fi) {
                    let gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[0]]) * this.tradeoff;
                    max_range = (360.0 - all_hues[arg_hues[fi]] + all_hues[arg_hues[0]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[0]]));
                }
                else {
                    let gap_hue = max(sqrt_nums[arg_hues[i]], sqrt_nums[arg_hues[i + 1]]) * this.tradeoff;
                    max_range = (all_hues[arg_hues[i + 1]] - all_hues[arg_hues[i]]) * (sqrt_nums[arg_hues[i]] / (gap_hue + sqrt_nums[arg_hues[i]] + sqrt_nums[arg_hues[i + 1]]));
                }
                let eq_range = min(min_range, max_range);
                if (modify_pcolors)
                    eq_range = max(eq_range, this.init_hue_min_delta);
                else
                    eq_range = max(eq_range, this.hue_min_range);

                if (!modify_pcolors && this.clusters[arg_hues[i]].num == 1) {
                    eq_range = 0.0;
                }

                this.clusters[arg_hues[i]].hue_min = fmod(all_hues[arg_hues[i]] - eq_range, 360.0);
                if (this.clusters[arg_hues[i]].hue_min < 0.0)
                    this.clusters[arg_hues[i]].hue_min += 360.0;

                this.clusters[arg_hues[i]].hue_max = fmod(all_hues[arg_hues[i]] + eq_range, 360.0);
                if (this.clusters[arg_hues[i]].hue_max < 0.0)
                    this.clusters[arg_hues[i]].hue_max += 360.0;

                this.clusters[arg_hues[i]].hue_range = eq_range;
            }
        }
        else {
            this.clusters[0].hue_min = 0.0;
            this.clusters[0].hue_max = 0.0;
            this.clusters[0].hue_range = 180;
            if (!modify_pcolors && !this.gen_color) {
                this.clusters[0].hue_range = min(this.estimate_delta * sqrt_nums[0], this.hue_max_range);
            }
        }

        let dist_matrix = basic_colors.map(() => basic_colors.map(() => 0));
        let min_basic_dist = 30.0;
        for (let i = 0; i < basic_colors.length; i++) {
            for (let j = i + 1; j < basic_colors.length; j++) {
                let dist = cdist(basic_colors[i], basic_colors[j]);
                dist_matrix[i][j] = dist;
                dist_matrix[j][i] = dist;
                min_basic_dist = min(min_basic_dist, dist);
            }
        }

        this.clusters.forEach((cluster, i) => {
            let min_dist = 10000.0;
            this.clusters.forEach((cluster2, j) => {
                if (i == j) return;
                let gap = max(sqrt_nums[i], sqrt_nums[j]) * this.tradeoff;
                let r = dist_matrix[i][j] * sqrt_nums[i] / (gap + sqrt_nums[i] + sqrt_nums[j]);
                if (r < min_dist)
                    min_dist = r;
            });
            if (min_dist > this.color_max_range) {
                if (!modify_pcolors && !this.gen_color) {
                    min_dist = min(this.estimate_delta * sqrt_nums[i], this.color_max_range);
                }
                else {
                    min_dist = this.color_max_range;
                }
            }
            if (modify_pcolors)
                min_dist = max(min_dist, this.init_color_min_delta);
            else
                min_dist = max(min_dist, this.color_min_range);
            if (!modify_pcolors && this.clusters[i].num == 1) {
                min_dist = 0.0;
            }
            this.clusters[i].range = min_dist;
        });

        // 2. discrimination: initialization and optimization
        this.shift_p = 1;
        let palette = [];
        this.clusters.forEach(cluster => {
            let range = [cluster.hue_range, 60, 60];
            if (cluster.num == 1) {
                palette.push(d3.rgb(cluster.color));
                return;
            }

            d3.range(cluster.num).forEach(() => {
                let color = d3.rgb(cluster.color);
                color = this.scope.disturbColor(color, range);
                if (modify_pcolors) {
                    color = modifyColorInConstraint(color, cluster.index, basic_colors);
                }
                color = modifyColorInColorRange(color, cluster.color, cluster.range, cluster.hue_min, cluster.hue_max);
                palette.push(color);
            });
        });

        let cur_score = this.evaluate_discrimination(palette)[0];
        let start_palette = palette.map(d => d);
        let start_score = cur_score;

        // trick: use blue noise to initialize and get estimated dart throw radius for normalization
        let dart_throw_r = 30.0;
        let lr = 0.25;
        if (basic_colors.length > 1) {
            dart_throw_r = max(min(min_basic_dist, dart_throw_r), this.global_color_dist);
        }
        let max_min_dist = 0.0;
        let ct = 0;
        while (dart_throw_r > 1e-5) {
            let find_satisfy = true;
            for (let i = 0; i < this.iter_num; i++) {
                let tp_palette = start_palette.map(d => d);
                let ct2 = 0, sign = 0;
                while ((sign = this.isDiscirminative(tp_palette, dart_throw_r)) != -1 && ct2 < 100) {
                    ct2 += 1;
                    let rcolor = this.scope.disturbColor(d3.rgb(tp_palette[sign]), [60, 60, 60]);
                    let cluster = this.clusters[this.idx2id[sign]];
                    if (modify_pcolors) {
                        rcolor = modifyColorInConstraint(rcolor, cluster.index, basic_colors);
                    }
                    rcolor = modifyColorInColorRange(rcolor, cluster.color, cluster.range, cluster.hue_min, cluster.hue_max);
                    tp_palette[sign] = rcolor;
                }
                if (sign != -1) {
                    find_satisfy = false;
                }
                let res = this.evaluate_discrimination(tp_palette);
                let new_score = res[0];
                let new_min_dist = res[1];
                if (new_score > cur_score) {
                    start_score = new_score;
                    start_palette = tp_palette;
                }
                if (new_min_dist > max_min_dist) {
                    max_min_dist = new_min_dist;
                }
            }
            if (find_satisfy || ct > 50) {
                break;
            }
            dart_throw_r -= (dart_throw_r - max_min_dist) * lr;
            ct += 1;
        }
        if (dart_throw_r > this.global_color_dist) {
            dart_throw_r = this.global_color_dist;
        }

        let best_palette = start_palette;
        let best_score = start_score;
        palette = start_palette;
        cur_score = start_score;
        let cur_temper = 1e-3;
        let end_temper = 1e-9;
        let dec = 0.99;
        let use_early_stop = true;
        let stop_iter_step0 = 50;
        let stop_iter = 10;
        let cur_stop_iter = 0;

        // to estimated up_bound and low_bound
        let harmony_score1 = -1 * template_score(palette);
        let harmony_score2 = -1 * geo_lc_linear(palette);
        let data_score = 0.0;
        if (use_data) data_score = this.evaluate_data(palette);
        let up_bound = [best_score, 0.0];
        let low_bound = [best_score, harmony_score1];
        let lc_bound = [harmony_score2, 0.0];
        let data_bound = [data_score, data_score];

        // optimization discrimination and estimate other bound
        const disturbPalette = (palette, disturb_step) => {
            if (randOne() < this.shift_p) {
                let new_palette = palette.map(d => d);
                let idx = randInt(0, new_palette.length - 1);
                let dis_color = this.scope.disturbColor(d3.rgb(new_palette[idx]), disturb_step);
                let cluster = this.clusters[this.idx2id[idx]];
                if (modify_pcolors) {
                    dis_color = modifyColorInConstraint(dis_color, cluster.index, basic_colors);
                }
                dis_color = modifyColorInColorRange(dis_color, cluster.color, cluster.range, cluster.hue_min, cluster.hue_max);
                new_palette[idx] = dis_color;
                return new_palette;
            }
            else {
                let new_palette = palette.map(d => d);
                let idx1 = randInt(0, new_palette.length - 1);
                let swap_part = this.idx2id[idx1];
                if (this.clusters[swap_part].num == 1) return new_palette;
                let idx2 = randInt(this.clusters[swap_part].start, this.clusters[swap_part].end - 1);
                let ct = 0;
                while (idx1 == idx2 && ct <= 5) {
                    idx2 = randInt(this.clusters[swap_part].start, this.clusters[swap_part].end - 1);
                    ct += 1;
                }
                let tmp = new_palette[idx1];
                new_palette[idx1] = new_palette[idx2];
                new_palette[idx2] = tmp;
                return new_palette;
            }

        };

        while (cur_temper > end_temper) {
            let tp_palette = palette;
            let tp_score = cur_score;
            let tp_best_score = best_score;
            for (let i = 0; i < this.iter_num; i++) {
                let new_palette = disturbPalette(palette, [30, 30, 30]);
                let res = this.evaluate_discrimination(new_palette);
                let new_score = res[0];
                let new_min_dist = res[1];
                let new_harmony_score1 = -1 * template_score(new_palette);
                let new_harmony_score2 = -1 * geo_lc_linear(new_palette);
                let new_data_score = 0.0;
                let delta_score = cur_score - new_score;
                if (use_data) this.evaluate_data(new_palette);

                if (new_score > best_score) {
                    best_score = new_score;
                    best_palette = new_palette;
                }
                if (delta_score <= 0 || randOne() < Math.exp((-delta_score) / cur_temper)) {
                    tp_palette = new_palette;
                    tp_score = new_score;
                }
                if (new_min_dist > dart_throw_r && new_score < low_bound[0]) {
                    low_bound[0] = new_score;
                }
                if (new_harmony_score1 < low_bound[1]) {
                    low_bound[1] = new_harmony_score1;
                }
                if (new_harmony_score2 < lc_bound[0]) {
                    lc_bound[0] = new_harmony_score2;
                }
                if (new_data_score < data_bound[0]) {
                    data_bound[0] = new_data_score;
                }
                if (new_data_score > data_bound[1]) {
                    data_bound[1] = new_data_score;
                }
            }
            palette = tp_palette;
            cur_score = tp_score;

            if (use_early_stop) {
                if (Math.abs(best_score - tp_best_score) < 1e-3) {
                    cur_stop_iter += 1;
                }
                else {
                    cur_stop_iter = 0;
                }
                if (cur_stop_iter >= stop_iter_step0) {
                    break;
                }
            }
            cur_temper *= dec;
        }
        let dist_score_without_penalty = this.evaluate_discrimination(best_palette, 0.1, 2.0, 0.0)[0];
        up_bound[0] = dist_score_without_penalty;
        up_bound[1] = 0.0;
        low_bound[0] = low_bound[0] - 1e-5;
        let tp_bound = [low_bound[1], up_bound[1]];
        this.harmony_args = [lc_bound[1] - lc_bound[0], this.bharg * (tp_bound[1] - tp_bound[0])];
        low_bound[1] = tp_bound[0] * this.harmony_args[0] + lc_bound[0] * this.harmony_args[1] - 1e-5;

        // 3. harmony: optimization harmony by dynamic weights
        cur_temper = 1e-3;
        end_temper = 1e-6;
        palette = best_palette;
        cur_score = best_score;
        let best_dis_score = this.evaluate_discrimination(palette)[0];
        let best_harmony_score = this.evaluate_harmony(palette);
        let start_dis_score = best_dis_score;
        let cur_score1 = best_dis_score;
        let cur_score2 = best_harmony_score;

        let batch_num = 5;
        let min_value = 0.1;
        let args = [1.0, 0.0];
        cur_stop_iter = 0;
        let cur_iter = 0;
        while (cur_temper > end_temper) {
            if (cur_iter % batch_num == 0) {
                let norm_score1 = (best_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
                let norm_score2 = (best_harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);
                let arg0 = 1.0 - norm_score1;
                let arg1 = 1.0 - norm_score2;
                if (arg0 < min_value) {
                    arg0 = min_value;
                }
                if (arg1 < min_value) {
                    arg1 = min_value;
                }
                let sum = arg0 + arg1;
                args[0] = arg0 / sum;
                args[1] = arg1 / sum;

                palette = best_palette;
                best_score = best_dis_score * args[0] + best_harmony_score * args[1];
                cur_score = best_score;
            }

            let tp_palette = palette;
            let tp_score = cur_score;
            let tp_best_score = best_score;
            let tp_dis_score = cur_score1;
            let tp_harmony_score = cur_score2;
            let norm_start_dis_score = (start_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
            let norm_tp_dis_score = (tp_dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
            let norm_tp_harmony_score = (tp_harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);

            for (let i = 0; i < this.iter_num; i++) {
                let new_palette = disturbPalette(palette, [10, 10, 10]);
                let dis_score = this.evaluate_discrimination(new_palette)[0];
                let harmony_score = this.evaluate_harmony(new_palette);
                let norm_dis_score = (dis_score - low_bound[0]) / (up_bound[0] - low_bound[0]);
                let norm_harmony_score = (harmony_score - low_bound[1]) / (up_bound[1] - low_bound[1]);
                let new_score = norm_dis_score * args[0] + norm_harmony_score * args[1];
                let delta_score = cur_score - new_score;
                let importance_factor = (norm_tp_dis_score > norm_tp_harmony_score) ?
                    (norm_dis_score > norm_harmony_score) : (norm_dis_score >= norm_start_dis_score);

                if (importance_factor) {
                    if (new_score > best_score) {
                        best_score = new_score;
                        best_palette = new_palette;
                        best_dis_score = dis_score;
                        best_harmony_score = harmony_score;
                    }
                    if (delta_score <= 0 || randOne() < Math.exp((-delta_score) / cur_temper)) {
                        tp_palette = new_palette;
                        tp_score = new_score;
                        tp_dis_score = dis_score;
                        tp_harmony_score = harmony_score;
                    }
                }
                if (dis_score > up_bound[0]) {
                    up_bound[0] = dis_score;
                }
            }
            palette = tp_palette;
            cur_score = tp_score;
            cur_score1 = tp_dis_score;
            cur_score2 = tp_harmony_score;

            if (use_early_stop) {
                if (Math.abs(tp_best_score - best_score) < 1e-3) {
                    cur_stop_iter += 1;
                }
                else {
                    cur_stop_iter = 0;
                }
                if (cur_stop_iter >= stop_iter) {
                    break;
                }
            }

            cur_temper *= dec;
            cur_iter += 1;
        }
        if (!use_data) return best_palette;

        // 4. data fitting: optimization data fitting by dynamic weights
        let alpha = args[1] / args[0];
        if (Math.abs(data_bound[1] - data_bound[0]) < 1e-5)
            data_bound[0] -= 1e-5;

        start_dis_score = best_dis_score;
        let start_harmony_score = best_harmony_score;
        let best_score_1 = start_dis_score + alpha * start_harmony_score;
        let best_score_2 = this.evaluate_data(best_palette);
        let dis_bound = [low_bound[0], up_bound[0]];
        let harmony_bound = [low_bound[1], up_bound[1]];

        cur_stop_iter = 0;
        cur_iter = 0;
        this.shift_p = 0.1;
        args[0] = 1.0;
        args[1] = 0.0;
        up_bound[0] = best_score_1;
        low_bound[0] += low_bound[1] * alpha;
        up_bound[1] = best_score_2;
        low_bound[1] = data_bound[0];

        cur_temper = 1e-3;
        end_temper = 1e-6;
        palette = best_palette;

        while (cur_temper > end_temper) {
            if (cur_iter % batch_num == 0) {
                let norm_score1 = (best_score_1 - low_bound[0]) / (up_bound[0] - low_bound[0]);
                let norm_score2 = (best_score_2 - low_bound[1]) / (up_bound[1] - low_bound[1]);
                let arg0 = 1.0 - norm_score1;
                let arg1 = 1.0 - norm_score2;
                if (arg0 < min_value) {
                    arg0 = min_value;
                }
                if (arg1 < min_value) {
                    arg1 = min_value;
                }
                let sum = arg0 + arg1;
                args[0] = arg0 / sum;
                args[1] = arg1 / sum;

                palette = best_palette;
                best_score = best_score_1 * args[0] + best_score_2 * args[1];
                cur_score = best_score;
            }
            let tp_palette = palette;
            let tp_score = cur_score;
            let tp_dis_score = this.evaluate_discrimination(palette)[0];
            let tp_harmony_score = this.evaluate_harmony(palette);
            let tp_data_score = this.evaluate_data(palette);

            let norm_start_dis_score = (start_dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
            let norm_start_harmony_score = (start_harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);
            let norm_tp_dis_score = (tp_dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
            let norm_tp_harmony_score = (tp_harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);
            let norm_tp_data_score = (tp_data_score - low_bound[1]) / (up_bound[1] - low_bound[1]);

            for (let i = 0; i < this.iter_num; i++) {
                let new_palette = disturbPalette(palette, [2, 2, 2]);
                let dis_score = this.evaluate_discrimination(new_palette)[0];
                let harmony_score = this.evaluate_harmony(new_palette);
                let score1 = dis_score + alpha * harmony_score;
                let score2 = this.evaluate_data(new_palette);

                let norm_score1 = (score1 - low_bound[0]) / (up_bound[0] - low_bound[0]);
                let norm_score2 = (score2 - low_bound[1]) / (up_bound[1] - low_bound[1]);
                let norm_dist = (dis_score - dis_bound[0]) / (dis_bound[1] - dis_bound[0]);
                let norm_harmony = (harmony_score - harmony_bound[0]) / (harmony_bound[1] - harmony_bound[0]);

                let new_score = norm_score1 * args[0] + norm_score2 * args[1];
                let delta_score = cur_score - new_score;
                let importance_factor = (norm_tp_dis_score > norm_tp_harmony_score) ?
                    (norm_dist > norm_harmony) || ((norm_tp_harmony_score > norm_tp_data_score) ? (norm_harmony > norm_score2) : (norm_harmony >= norm_start_harmony_score))
                    : (norm_dist >= norm_start_dis_score);

                if (importance_factor) {
                    if (new_score > best_score) {
                        best_score = new_score;
                        best_palette = new_palette;
                        best_score_1 = score1;
                        best_score_2 = score2;
                    }
                    if (delta_score <= 0 || randOne() < Math.exp((-delta_score) / cur_temper)) {
                        tp_palette = new_palette;
                        tp_score = new_score;
                    }
                }
                if (score2 > up_bound[1]) {
                    up_bound[1] = score2;
                }
            }
            palette = tp_palette;
            cur_score = tp_score;

            cur_iter += 1;
            cur_temper *= dec;
        }
        return best_palette;
    }
    isDiscirminative(palette, r) {
        let idx = -1;
        let n = palette.length;
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                let dist = cdist(palette[i], palette[j]);
                if (dist < r) {
                    idx = i;
                    break;
                }
            }
            if (idx != -1) break;
        }
        return idx;
    }

    evaluate_discrimination(palette, lambda1 = 0.1, lambda2 = 2.0, lambda3 = 1.0) {
        // 1. min color distance
        let n = palette.length;
        let min_dist = 1000.0;
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                let dist = cdist(palette[i], palette[j]);
                min_dist = min(min_dist, dist);
            }
        }
        // 2. avg name difference
        let avg_name_diff = 0.0;
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                avg_name_diff += getNameDifference(palette[i], palette[j]);
            }
        }
        avg_name_diff /= n * (n - 1) / 2;

        // 3. penalty
        let penalty = min(min_dist - 10, 0);
        return [min_dist * lambda1 + avg_name_diff * lambda2 + penalty * lambda3, min_dist];
    }

    evaluate_harmony(palette) {
        // 1. template score
        let score = template_score(palette);
        // 2. geometric lc
        let geo_lc = geo_lc_linear(palette);
        return -score * this.harmony_args[0] - geo_lc * this.harmony_args[1];
    }

    load_data(data) {
        this.data = data;
        // TODO: implement this function for specific data format
    }

    evaluate_data(palette) {
        return 0.0;
        // TODO: implement this function for specific data format and mode(such as discriminative mode, harmony mode, etc.)
    }
}

function main() {
    // use case for testing
    let dc = new DynamicColor();
    let basic_colors = [[255, 0, 0], [0, 255, 0], [0, 0, 255]];
    let palette_sizes = [5, 5, 5];
    let modify_pcolors = false;
    let res = dc.run([], palette_sizes, true);
    debugger;
}