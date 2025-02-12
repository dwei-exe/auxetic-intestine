const Eng = {};

Eng.Vector = class {
  constructor(vals) {
    this.v = vals;
  }
  add(vec) {
    for (let i = 0; i < this.v.length; i++) {
      this.v[i] += vec[i];
    }
  }
  sub(vec) {
    for (let i = 0; i < this.v.length; i++) {
      this.v[i] -= vec[i];
    }
  }
  mult(n) {
    for (let i = 0; i < this.v.length; i++) {
      this.v[i] *= n;
    }
  }
  div(n) {
    for (let i = 0; i < this.v.length; i++) {
      this.v[i] /= n;
    }
  }
  get x() {
    return this.v[0];
  }
  get y() {
    return this.v[1];
  }
  get z() {
    return this.v[2];
  }
  set x(new_x) {
    this.v[0] = new_x;
  }
  set y(new_y) {
    this.v[1] = new_y;
  }
  set z(new_z) {
    this.v[2] = new_z;
  }
  get_spherical() {
    let theta_1 = Math.atan2(this.z, this.x);
    let theta_2 = Math.atan2(this.y, Math.sqrt(this.z ** 2 + this.y ** 2));
    return new Eng.Vector([theta_1, theta_2]);
  }
  get_Euclidean_length() {
    return Math.sqrt(Eng.Vector.dot(this.v, this.v));
  }
  get_clone() {
    return new Eng.Vector([...this.v]);
  }
};
Eng.Vector.dot = function (vec1, vec2) {
  let sum = 0;
  for (let i = 0; i < vec1.length; i++) {
    sum += vec1[i] * vec2[i];
  }
  return sum;
};
Eng.Vector.normalized = function (vec) {
  let sum = 0;
  for (let i = 0; i < vec.length; i++) {
    sum += vec[i] ** 2;
  }
  let mag = Math.sqrt(sum);

  let normalized = [];
  for (let i = 0; i < vec.length; i++) {
    normalized[i] = vec[i] / mag;
  }
  return new Eng.Vector(normalized);
};
Eng.Vector.average = function (vecList, vNum) {
  let sum = 0;
  for (let i = 0; i < vecList.length; i++) {
    sum += vecList[i].v[vNum];
  }
  sum /= vecList.length;
  return sum;
};
Eng.Vector.add = function (vec1, vec2) {
  let vec3 = new Eng.Vector([]);
  for (let i = 0; i < vec1.length; i++) {
    vec3.v[i] = vec1[i] + vec2[i];
  }
  return vec3;
};
Eng.Vector.sub = function (vec1, vec2) {
  let vec3 = new Eng.Vector([]);
  for (let i = 0; i < vec1.length; i++) {
    vec3.v[i] = vec1[i] - vec2[i];
  }
  return vec3;
};
Eng.Vector.mult = function (vec1, s1) {
  let vec3 = new Eng.Vector([]);
  for (let i = 0; i < vec1.length; i++) {
    vec3.v[i] = vec1[i] * s1;
  }
  return vec3;
};
Eng.Vector.div = function (vec1, s1) {
  let vec3 = new Eng.Vector([]);
  for (let i = 0; i < vec1.length; i++) {
    vec3.v[i] = vec1[i] / s1;
  }
  return vec3;
};
Eng.Vector.cross = function (vec1, vec2) {
  let vec3 = new Eng.Vector([]);
  vec3.v = [vec1[1] * vec2[2] - vec1[2] * vec2[1], vec1[2] * vec2[0] - vec1[0] * vec2[2], vec1[0] * vec2[1] - vec1[1] * vec2[0]];
  return vec3;
};

Eng.Matrix = class {
  constructor(rows, cols) {
    this.m = [];
    for (let r = 0; r < rows; r++) {
      this.m[r] = [];
      for (let c = 0; c < cols; c++) {
        this.m[r][c] = 0;
      }
    }
  }
};
Eng.Matrix.multMatrix = function (mat1, mat2) {
  if (mat1[0].length !== mat2.length) return false;

  let mat3 = new Eng.Matrix(mat1.length, mat2[0].length);

  for (let r = 0; r < mat1.length; r++) {
    for (let c = 0; c < mat2[0].length; c++) {
      for (let i = 0; i < mat1[0].length; i++) {
        mat3.m[r][c] += mat1[r][i] * mat2[i][c];
      }
    }
  }

  return mat3;
};
Eng.Matrix.multScalar = function (mat1, s1) {
  let mat2 = new Eng.Matrix(mat1.length, mat1[0].length);

  for (let r = 0; r < mat1.length; r++) {
    for (let c = 0; c < mat1[0].length; c++) {
      mat2.m[r][c] = mat1[r][c] * s1;
    }
  }

  return mat2;
};
Eng.Matrix.addMatrix = function (mat1, mat2) {
  if (mat1.length !== mat2.length || mat1[0].length !== mat2[0].length) return false;

  let mat3 = new Eng.Matrix(mat1.length, mat1[0].length);

  for (let r = 0; r < mat1.length; r++) {
    for (let c = 0; c < mat1[0].length; c++) {
      mat3.m[r][c] = mat1[r][c] + mat2[r][c];
    }
  }

  return mat3;
};
Eng.Matrix.vectToMat = function (vec) {
  let mat = new Eng.Matrix(vec.length, 1);

  for (let i = 0; i < vec.length; i++) {
    mat.m[i][0] = vec[i];
  }

  return mat;
};
Eng.Matrix.vectToMat2 = function (vec) {
  let mat = new Eng.Matrix(1, vec.length);

  for (let i = 0; i < vec.length; i++) {
    mat.m[0][i] = vec[i];
  }

  return mat;
};
Eng.Matrix.matToVect = function (mat) {
  let vec = new Eng.Vector([]);

  for (let i = 0; i < mat.length; i++) {
    vec.v[i] = mat[i][0];
  }

  return vec;
};
Eng.Matrix.projection = function (cam) {
  let aspectRatio = cam.a;
  let theta = cam.theta;
  let far = cam.far;
  let near = cam.near;

  let fov = 1 / Math.tan(theta / 2);
  let zNorm = far / (far - near);

  let projectionMatrix = new Eng.Matrix(4, 4);
  projectionMatrix.m[0][0] = aspectRatio * fov;
  projectionMatrix.m[1][1] = fov;
  projectionMatrix.m[2][2] = zNorm;
  projectionMatrix.m[3][2] = -near * zNorm;
  projectionMatrix.m[2][3] = 1;

  return projectionMatrix;
};
Eng.Matrix.pointAt = function (cam, upE) {
  let camPos = cam.pos;
  let lookDir = cam.rot;

  let target = new Eng.Vector([0, 0, 1]);
  let rotationYMatrix = Eng.Matrix.rotationY(lookDir.v);
  target = Eng.Matrix.vectToMat(target.v);
  target = Eng.Matrix.multMatrix(rotationYMatrix.m, target.m);
  target = Eng.Matrix.matToVect(target.m);
  target.add(camPos.v);

  let forward = Eng.Vector.sub(target.v, camPos.v);
  let normForward = Eng.Vector.normalized(forward.v);

  let dot = Eng.Vector.dot(upE.v, normForward.v);
  let a = Eng.Vector.mult(normForward.v, dot);
  let up = Eng.Vector.sub(upE.v, a.v);
  let normUp = Eng.Vector.normalized(up.v);

  let normRight = Eng.Vector.cross(normUp.v, normForward.v);

  let pointAtMatrix = new Eng.Matrix(4, 4);
  pointAtMatrix.m[0][0] = normRight.v[0];
  pointAtMatrix.m[0][1] = normRight.v[1];
  pointAtMatrix.m[0][2] = normRight.v[2];
  pointAtMatrix.m[1][0] = normUp.v[0];
  pointAtMatrix.m[1][1] = normUp.v[1];
  pointAtMatrix.m[1][2] = normUp.v[2];
  pointAtMatrix.m[2][0] = normForward.v[0];
  pointAtMatrix.m[2][1] = normForward.v[1];
  pointAtMatrix.m[2][2] = normForward.v[2];
  pointAtMatrix.m[3][0] = camPos.v[0];
  pointAtMatrix.m[3][1] = camPos.v[1];
  pointAtMatrix.m[3][2] = camPos.v[2];
  pointAtMatrix.m[3][3] = 1;

  return pointAtMatrix;
};
Eng.Matrix.lookAt = function (pointAtMatrix) {
  let pointInverse = new Eng.Matrix(4, 4);
  pointInverse.m[0][0] = pointAtMatrix[0][0];
  pointInverse.m[0][1] = pointAtMatrix[1][0];
  pointInverse.m[0][2] = pointAtMatrix[2][0];
  pointInverse.m[1][0] = pointAtMatrix[0][1];
  pointInverse.m[1][1] = pointAtMatrix[1][1];
  pointInverse.m[1][2] = pointAtMatrix[2][1];
  pointInverse.m[2][0] = pointAtMatrix[0][2];
  pointInverse.m[2][1] = pointAtMatrix[1][2];
  pointInverse.m[2][2] = pointAtMatrix[2][2];
  pointInverse.m[3][0] = -(
    pointAtMatrix[3][0] * pointInverse.m[0][0] +
    pointAtMatrix[3][1] * pointInverse.m[1][0] +
    pointAtMatrix[3][2] * pointInverse.m[2][0]
  );
  pointInverse.m[3][1] = -(
    pointAtMatrix[3][0] * pointInverse.m[0][1] +
    pointAtMatrix[3][1] * pointInverse.m[1][1] +
    pointAtMatrix[3][2] * pointInverse.m[2][1]
  );
  pointInverse.m[3][2] = -(
    pointAtMatrix[3][0] * pointInverse.m[0][2] +
    pointAtMatrix[3][1] * pointInverse.m[1][2] +
    pointAtMatrix[3][2] * pointInverse.m[2][2]
  );
  pointInverse.m[3][3] = 1;

  return pointInverse;
};
Eng.Matrix.rotationX = function (v) {
  let rotationXMatrix = new Eng.Matrix(3, 3);
  rotationXMatrix.m[0][0] = 1;
  rotationXMatrix.m[1][1] = Math.cos(v[0]);
  rotationXMatrix.m[1][2] = -Math.sin(v[0]);
  rotationXMatrix.m[2][1] = Math.sin(v[0]);
  rotationXMatrix.m[2][2] = Math.cos(v[0]);

  return rotationXMatrix;
};
Eng.Matrix.rotationY = function (v) {
  let rotationYMatrix = new Eng.Matrix(3, 3);
  rotationYMatrix.m[0][0] = Math.cos(v[1]);
  rotationYMatrix.m[0][2] = Math.sin(v[1]);
  rotationYMatrix.m[1][1] = 1;
  rotationYMatrix.m[2][0] = -Math.sin(v[1]);
  rotationYMatrix.m[2][2] = Math.cos(v[1]);

  return rotationYMatrix;
};
Eng.Matrix.rotationZ = function (v) {
  let rotationZMatrix = new Eng.Matrix(3, 3);
  rotationZMatrix.m[0][0] = Math.cos(v[2]);
  rotationZMatrix.m[0][1] = -Math.sin(v[2]);
  rotationZMatrix.m[1][0] = Math.sin(v[2]);
  rotationZMatrix.m[1][1] = Math.cos(v[2]);
  rotationZMatrix.m[2][2] = 1;

  return rotationZMatrix;
};

Eng.Scene = class {
  constructor() {
    this.meshes = [];
    this.lights = [];
  }
  addMesh(mesh) {
    this.meshes.push(mesh);
  }
  addLight(light) {
    this.lights.push(light);
  }
  findById(id) {
    for (let mesh of this.meshes) {
      if (mesh.id === id) {
        return mesh;
      }
    }
    for (let light of this.lights) {
      if (light.id === id) {
        return light;
      }
    }
  }
};

Eng.Camera = class {
  constructor(pos, rot, a, theta, far, near) {
    this.pos = pos;
    this.rot = rot;

    this.a = a;
    this.theta = theta;
    this.far = far;
    this.near = near;
  }
  translate(vec) {
    this.pos.add(vec.v);
  }
  rotate(vec) {
    this.rot.add(vec.v);
  }
};

Eng.DirectionalLight = class {
  constructor(rot) {
    this.rot = rot;
  }
};

Eng.Polygon = class {
  constructor(vecs) {
    this.vecs = vecs;
  }
};

Eng.Mesh = class {
  constructor(pols, id) {
    this.pols = pols;
    this.id = id;
  }
  translate(v) {
    for (let i = 0; i < this.pols.length; i++) {
      let pol = this.pols[i].vecs;
      for (let j = 0; j < pol.length; j++) {
        this.pols[i].vecs[j].add(v.v);
      }
    }
  }
  rotate(v) {
    let rotationXMatrix = Eng.Matrix.rotationX(v.v);
    let rotationYMatrix = Eng.Matrix.rotationY(v.v);
    let rotationZMatrix = Eng.Matrix.rotationZ(v.v);

    for (let i = 0; i < this.pols.length; i++) {
      for (let j = 0; j < this.pols[i].vecs.length; j++) {
        let vec = this.pols[i].vecs[j];

        let vectMatrix = Eng.Matrix.vectToMat(vec.v);

        vectMatrix = Eng.Matrix.multMatrix(rotationXMatrix.m, vectMatrix.m);
        vectMatrix = Eng.Matrix.multMatrix(rotationYMatrix.m, vectMatrix.m);
        vectMatrix = Eng.Matrix.multMatrix(rotationZMatrix.m, vectMatrix.m);

        this.pols[i].vecs[j] = Eng.Matrix.matToVect(vectMatrix.m);
      }
    }
  }
};

Eng.Renderer = class {
  constructor() {
    this.canvas = document.createElement("canvas");
    this.canvas.width = 500;
    this.canvas.height = 500;

    this.ctx = this.canvas.getContext("2d");

    this.backgroundColor = "grey";
    this.faceColor = "white";
    this.edgeColor = "black";
  }
  appendCanvas(element) {
    element.appendChild(this.canvas);
  }
  render(scene, cam) {
    let pointAtMatrix = Eng.Matrix.pointAt(cam, new Eng.Vector([0, 1, 0]));

    let lookAtMatrix = Eng.Matrix.lookAt(pointAtMatrix.m);

    let projectionMatrix = Eng.Matrix.projection(cam);

    let meshes = scene.meshes;

    let polygons = [];

    for (let mesh of meshes) {
      for (let pol of mesh.pols) {
        polygons.push(new Eng.Polygon(pol.vecs));
      }
    }

    polygons.sort(function (a, b) {
      return Eng.Vector.average(a.vecs, 2) - Eng.Vector.average(b.vecs, 2);
    });

    let lights = scene.lights;

    let projectedPolygons = [];

    for (let pol of polygons) {
      let lines = [];
      lines[0] = Eng.Vector.sub(pol.vecs[1].v, pol.vecs[0].v);
      lines[1] = Eng.Vector.sub(pol.vecs[2].v, pol.vecs[0].v);

      let crossVector = Eng.Vector.cross(lines[0].v, lines[1].v);

      let normalizedCross = Eng.Vector.normalized(crossVector.v);

      let projectedVectors = [];

      for (let vec of pol.vecs) {
        let vectMat = Eng.Matrix.vectToMat2(vec.v);
        vectMat.m[0][3] = 1;
        let transMatrix = Eng.Matrix.multMatrix(vectMat.m, lookAtMatrix.m);
        transMatrix = Eng.Matrix.vectToMat(transMatrix.m[0]);

        let projectedVectorsMatrix = Eng.Matrix.multMatrix(projectionMatrix.m, transMatrix.m);

        let dividedXYZMatrix = Eng.Matrix.multScalar(projectedVectorsMatrix.m, 1 / projectedVectorsMatrix.m[3][0]);

        let addMatrix = new Eng.Matrix(4, 1);
        addMatrix.m[0][0] = 1;
        addMatrix.m[1][0] = 1;

        let addedOneMatrix = Eng.Matrix.addMatrix(dividedXYZMatrix.m, addMatrix.m);

        let scaleMatrix = new Eng.Matrix(4, 4);
        scaleMatrix.m[0][0] = this.canvas.width / 2;
        scaleMatrix.m[1][1] = this.canvas.height / 2;
        scaleMatrix.m[2][2] = 1;
        scaleMatrix.m[3][3] = 1;

        let scaledVectorsMatrix = Eng.Matrix.multMatrix(scaleMatrix.m, addedOneMatrix.m);

        projectedVectors.push(Eng.Matrix.matToVect(scaledVectorsMatrix.m));
      }

      let flatPolygon = new Eng.Polygon(projectedVectors);

      let polLight = 0;
      for (let light of lights) {
        polLight += Eng.Vector.dot(light.rot.v, normalizedCross.v);
      }

      flatPolygon.color = `hsl(200, 100%, ${-polLight * 100}%)`;

      projectedPolygons.push(flatPolygon);
    }

    this.ctx.fillStyle = this.backgroundColor;
    this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

    for (let pol of projectedPolygons) {
      this.ctx.strokeStyle = pol.color;
      this.ctx.fillStyle = pol.color;
      this.ctx.beginPath();
      this.ctx.moveTo(pol.vecs[0].v[0], pol.vecs[0].v[1]);
      for (let i = 1; i < pol.vecs.length; i++) {
        this.ctx.lineTo(pol.vecs[i].v[0], pol.vecs[i].v[1]);
      }
      this.ctx.closePath();
      this.ctx.stroke();
      this.ctx.fill();
    }

    this.ctx.fillStyle = "red";
    this.ctx.font = "50px Arial";
    this.ctx.fillText("DINOSAR", 100, 50);
  }
};

Eng.Readline = function (text) {
  let localString = "";
  for (let i = Eng.Readline.index; i < text.length; i++) {
    localString += text[i];
    if (text[i] == "\n") {
      i++;
      Eng.Readline.index = i;
      return localString;
    }
  }
  return false;
};
Eng.Readline.index = 0;

Eng.ObjLoader = function (text) {
  let vertices = [];

  let pols = [];

  let name;

  while (true) {
    let line = Eng.Readline(text);

    if (line === false) break;

    if (line[0] === "#") continue;

    if (line[0] === "o") name = line.slice(1);

    if (line[0] === "v") {
      let coords = [];

      let localNum = "";

      for (let i = 2; i < line.length; i++) {
        if (line[i] === " " || line[i] == "\n") {
          coords.push(localNum);
          localNum = "";
        } else {
          localNum += line[i];
        }
      }

      let vertex = [];
      vertex[0] = coords[0];
      vertex[1] = coords[1];
      vertex[2] = coords[2];

      vertices.push(vertex);
    }

    if (line[0] === "s") continue;

    if (line[0] === "f") {
      let pol = [];

      let localNum = "";

      for (let i = 2; i < line.length; i++) {
        if (line[i] === " " || line[i] == "\n") {
          pol.push(vertices[parseInt(localNum) - 1]);
          localNum = "";
        } else {
          localNum += line[i];
        }
      }

      pols.push(pol);
    }
  }

  return pols;
};

Eng.ParseMesh = function (text) {
  let mesh = Eng.ObjLoader(text);
  let pols = [];
  for (let pol of mesh) {
    let vectors = [];
    for (let vector of pol) {
      let vectorE = new Eng.Vector([parseFloat(vector[0]), parseFloat(vector[1]), parseFloat(vector[2])]);
      vectors.push(vectorE);
    }
    let polE = new Eng.Polygon(vectors);
    pols.push(polE);
  }
  let meshE = new Eng.Mesh(pols);

  return meshE;
};

export default Eng;
