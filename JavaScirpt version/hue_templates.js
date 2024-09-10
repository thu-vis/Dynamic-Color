// hue harmony template implementation
// matsuda_templates: Matsuda et al. 2012
// geometric_hue_templates: geometric hue harmony

// important tip: matsuda templates use hsl hue, while geometric templates use hcl hue
var _norm_h = function (value) {
    return (value % 360 + 360) % 360;
};

class hue_templates {
    constructor(ranges) {
        let that = this;
        ranges.forEach(range => {
            range[0] = _norm_h(range[0]);
            range[1] = _norm_h(range[1]);
        });
        that.ranges = ranges;
        that.key_positions = [];
        that.ranges.forEach(range => {
            that.key_positions.push(range[0]);
            that.key_positions.push(range[1]);
            that.key_positions.push(_norm_h((range[0] + range[1]) / 2));
        });
    };

    rotate(alpha) {
        return new hue_templates(this.ranges.map(range => {
            return [_norm_h(range[0] + alpha), _norm_h(range[1] + alpha)];
        }))
    };

    mindist(hue) {
        let mindist = undefined;
        for (let i = 0; i < this.ranges.length; i++) {
            let range = this.ranges[i];
            let dist;
            if (range[0] < range[1]) {
                if (hue < range[0]) {
                    dist = Math.min(range[0] - hue, hue + 360 - range[1]);
                }
                else if (hue > range[1]) {
                    dist = Math.min(hue - range[1], range[0] + 360 - hue);
                }
                else return 0;
            }
            if (range[1] <= range[0]) {
                if (hue >= range[0] || hue <= range[1]) {
                    return 0;
                }
                dist = Math.min(range[0] - hue, hue - range[1]);
            }
            if (mindist === undefined || dist < mindist) {
                mindist = dist;
            }
        };
        return mindist;
    }
};

const check_color_min_dist = function (template, hues) {
    let rotate_list = [];
    let rotate_set = new Set();
    hues.forEach(hue => {
        hue = _norm_h(hue);
        template.key_positions.forEach(key => {
            let r = _norm_h(hue - key);
            if (!rotate_set.has(r)) {
                rotate_list.push(r);
                rotate_set.add(r);
            }
        });
    });

    let min_dist = undefined;
    let min_dist_rotate = undefined;
    rotate_list.forEach(r => {
        let cur_template = template.rotate(r);
        let dist = 0;
        hues.forEach(hue => {
            dist += cur_template.mindist(hue);
        });
        if (min_dist === undefined || dist < min_dist) {
            min_dist = dist;
            min_dist_rotate = r;
        }
    });
    return [min_dist, min_dist_rotate];
}

class matsuda_templates {
    constructor() {
        let that = this;
        that.template_ranges = {
            'i': [[-9, 9]],
            'V': [[-46.8, 46.8]],
            'L': [[-9, 9], [50.4, 129.6]],
            'I': [[-9, 9], [171, 189]],
            'T': [[0, 180]],
            'Y': [[-46.8, 46.8], [171, 189]],
            'X': [[-46.8, 46.8], [133.2, 226.8]],
        };
        that.use_names = ['i', 'V', 'L', 'I', 'T', 'Y', 'X'];
        that.templates = that.use_names.map(name => {
            let ranges = that.template_ranges[name];
            return new hue_templates(ranges);
        });
    };

    check_templates(hues) {
        // should be hsv hues
        let min_dist = undefined;
        let min_dist_rotate = undefined;
        let min_dist_template_id = undefined;
        this.templates.forEach((template, id) => {
            let [dist, rotate] = check_color_min_dist(template, hues);
            if (min_dist === undefined || dist < min_dist) {
                min_dist = dist;
                min_dist_rotate = rotate;
                min_dist_template_id = id;
            }
        });
        return [min_dist, this.use_names[min_dist_template_id], this.templates[min_dist_template_id], min_dist_rotate];
    }


};

class geometric_hue_templates {
    constructor() {
        let that = this;
        that.template_ranges = {
            'analog': [[-45, 45]],
            'complimentary': [[-15, 15], [165, 195]],
            'triad': [[-15, 15], [105, 135], [225, 255]],
        };
        that.use_names = ['analog', 'complimentary', 'triad'];
        that.templates = that.use_names.map(name => {
            let ranges = that.template_ranges[name];
            return new hue_templates(ranges);
        });
    };

    check_templates(hues) {
        // should be hsv hues
        let min_dist = undefined;
        let min_dist_rotate = undefined;
        let min_dist_template_id = undefined;
        this.templates.forEach((template, id) => {
            let [dist, rotate] = check_color_min_dist(template, hues);
            if (min_dist === undefined || dist < min_dist) {
                min_dist = dist;
                min_dist_rotate = rotate;
                min_dist_template_id = id;
            }
        });
        return [min_dist, this.use_names[min_dist_template_id], this.templates[min_dist_template_id], min_dist_rotate];
    }


};

module.exports = matsuda_templates;
