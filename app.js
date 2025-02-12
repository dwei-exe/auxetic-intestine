import url from "url";
import path from "path";
import http from "http";

import express from "express";

import log from "./lib/logger.js";

const __dirname = path.dirname(url.fileURLToPath(import.meta.url));

const app = express();
app.use("/public", express.static(path.join(__dirname, "public")));

app.get("/mathjs", function (req, res) {
  res.sendFile(path.join(__dirname, "node_modules/mathjs/lib/browser/math.js"));
});
app.get("/", function (req, res) {
  res.sendFile(path.join(__dirname, "index.html"));
});

const server = http.createServer(app);
server.listen(process.env.PORT || 80, function () {
  log("Server", "listening");
});
