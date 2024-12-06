import fs from "fs";
import path from "path";
import url from "url";
import { createCanvas } from "canvas";

import * as LINALG from "./linear-calc.js";
import { Point_Dist } from "./point-dist.js";

const colon = new LINALG.Parametric2D(
  (s, t) => (s + Math.pi) / 3,
  (s, t) => Math.cos(t) * (1 + (1 / 5) * sin(s)),
  (s, t) => Math.sin(t) * (1 + (1 / 5) * sin(s)),
  [-20, 20],
  [0, 2 * Math.PI]
);

const resolution_d = 1;
const point_dist = new Point_Dist(colon, resolution_d);

{
  const print_scale = 10;
  const canvas = createCanvas(
    (colon.s_range[1] - colon.s_range[0]) * print_scale,
    (colon.t_range[1] - colon.t_range[0]) * print_scale
  );
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = "white";
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  for (let point of point_dist.points) {
    ctx.fillStyle = "#000000";
    ctx.beginPath();
    ctx.arc(
      (point.x - colon.s_range[0]) * print_scale,
      (point.y - colon.t_range[0]) * print_scale,
      (resolution_d * print_scale) / 4,
      0,
      2 * Math.PI
    );
    ctx.fill();
  }

  const __dirname = path.dirname(url.fileURLToPath(import.meta.url));

  const out = fs.createWriteStream(path.join(__dirname, "test.png"));
  const stream = canvas.createPNGStream();
  stream.pipe(out);
  out.on("finish", function () {
    console.log("PNG exported.");
  });
}
