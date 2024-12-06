import { Vector2D } from "./linear-calc.js";

const Point_Dist = class {
  constructor(parametric2D, resolution_d) {
    this.points = [];
    for (
      let s = parametric2D.s_range[0];
      s < parametric2D.s_range[1];
      s += resolution_d
    ) {
      for (
        let t = parametric2D.t_range[0];
        t < parametric2D.t_range[1];
        t += resolution_d
      ) {
        this.points.push(new Vector2D(s, t));
      }
    }
  }
};

export { Point_Dist };
