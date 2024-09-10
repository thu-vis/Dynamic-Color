var d3 = require("./lib/d3.v7.min.js");
var utils = require('./lib/utils.js');

// tips: 40 for original visualization, 50 for lighter environment
const colorscope = {
    'hue_scope': [0, 360],
    'chroma_scope': [50, 85],
    'lumi_scope': [50, 85],
    'all_hue': false,
};

// color scope to control the color range and disturb colors
class ColorScope {
    constructor(scope = colorscope) {
        this.scope = scope;
        this.randf = utils.randomDouble;
        this.margin = 5;
    }

    setMode(center = false) {
        if (center) {
            this.scope['lumi_scope'] = [55, 80];
            this.scope['chroma_scope'] = [55, 80];
        }
        else {
            // the same with the previous setting
            this.scope['lumi_scope'] = [50, 85];
            this.scope['chroma_scope'] = [50, 85];
        }
    }

    randomHue() {
        if (this.scope['all_hue']) {
            return this.randf(0, 360);
        }
        let hue = this.randf(this.scope['hue_scope'][0], this.scope['hue_scope'][1]);
        if (hue > 84 && hue < 115) {
            if (hue < 99.5) hue = 84;
            else hue = 115;
        }
        return hue;
    }

    disturbColor(color, randg) {
        // randf: random function, randg: random scope [3]
        let hcl = d3.hcl(color);
        let hue = hcl.h;
        let chroma = this.randf(Math.max(Math.round(hcl.c - randg[1]), this.scope['chroma_scope'][0] + this.margin),
            Math.min(Math.round(hcl.c + randg[1]), this.scope['chroma_scope'][1]));
        let lumi = this.randf(Math.max(Math.round(hcl.l - randg[2]), this.scope['lumi_scope'][0] + this.margin),
            Math.min(Math.round(hcl.l + randg[2]), this.scope['lumi_scope'][1]));

        if (this.scope['all_hue']) {
            hue = hcl.h + this.randf(-randg[0], randg[0]);
            hue = hue % 360 + (hue < 0 ? 360 : 0);
        } else {
            hue = this.randf(Math.max(Math.round(hcl.h - randg[0]), this.scope['hue_scope'][0]),
                Math.min(Math.round(hcl.h + randg[0]), this.scope['hue_scope'][1]));
        }
        if (hue > 84 && hue < 115) {
            if (hue < 99.5) hue = 84;
            else hue = 115;
        }

        return d3.hcl(hue, chroma, lumi);
    }

    modifyColorInScope(color) {
        let hcl = d3.hcl(color);
        let hue = hcl.h;
        let chroma = Math.max(this.scope['chroma_scope'][0], Math.min(hcl.c, this.scope['chroma_scope'][1]));
        let lumi = Math.max(this.scope['lumi_scope'][0], Math.min(hcl.l, this.scope['lumi_scope'][1]));

        if (this.scope['all_hue']) {
            hue = hue % 360 + (hue < 0 ? 360 : 0);
        } else {
            hue = Math.max(this.scope['hue_scope'][0], Math.min(hcl.h, this.scope['hue_scope'][1]));
        }
        if (hue > 84 && hue < 115) {
            if (hue < 99.5) hue = 84;
            else hue = 115;
        }

        return d3.hcl(hue, chroma, lumi);
    }
}

module.exports = ColorScope;
