const Vector3D = class {
  constructor(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  length() {
    return Math.sqrt(this.x ** 2 + this.y ** 2 + this.z ** 2);
  }
};
Vector3D.difference = function (v1, v2) {
  return new Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
};
Vector3D.cross = function (v1, v2) {
  return new Vector3D(
    v1.y * v2.z - v1.z * v2.y,
    v1.z * v2.x - v1.x * v2.z,
    v1.x * v2.y - v1.y * v2.x
  );
};

const Vector2D = class {
  constructor(x, y) {
    this.x = x;
    this.y = y;
  }
};

const Matrix = class {
  constructor(row, col, builder_func) {
    this.data = [];
    for (let i = 0; i < row; i++) {
      this.data[i] = [];
      for (let j = 0; j < col; j++) {
        this.data[i][j] = builder_func(i, j);
      }
    }
  }
};

const Parametric2D = class {
  constructor(x_func, y_func, z_func, s_range, t_range) {
    this.x_func = x_func;
    this.y_func = y_func;
    this.z_func = z_func;

    this.s_range = s_range;
    this.t_range = t_range;
  }
  expansion(s, t, d) {
    let dv_left = Vector3D.difference(
      new Vector3D(
        this.x_func(s - d, t),
        this.y_func(s - d, t),
        this.z_func(s - d, t)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let dv_right = Vector3D.difference(
      new Vector3D(
        this.x_func(s + d, t),
        this.y_func(s + d, t),
        this.z_func(s + d, t)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let dv_down = Vector3D.difference(
      new Vector3D(
        this.x_func(s, t - d),
        this.y_func(s, t - d),
        this.z_func(s, t - d)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let dv_up = Vector3D.difference(
      new Vector3D(
        this.x_func(s, t + d),
        this.y_func(s, t + d),
        this.z_func(s, t + d)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let upper_right = Vector3D.cross(dv_right, dv_up).length();
    let upper_left = Vector3D.cross(dv_up, dv_left).length();
    let lower_left = Vector3D.cross(dv_left, dv_down).length();
    let lower_right = Vector3D.cross(dv_down, dv_right).length();

    let expansion = upper_right + upper_left + lower_left + lower_right;

    return expansion;
  }
  curving(s, t, d) {
    let dv_left = Vector3D.difference(
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t)),
      new Vector3D(
        this.x_func(s - d, t),
        this.y_func(s - d, t),
        this.z_func(s - d, t)
      )
    );

    let dv_right = Vector3D.difference(
      new Vector3D(
        this.x_func(s + d, t),
        this.y_func(s + d, t),
        this.z_func(s + d, t)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let dv_down = Vector3D.difference(
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t)),
      new Vector3D(
        this.x_func(s, t - d),
        this.y_func(s, t - d),
        this.z_func(s, t - d)
      )
    );

    let dv_up = Vector3D.difference(
      new Vector3D(
        this.x_func(s, t + d),
        this.y_func(s, t + d),
        this.z_func(s, t + d)
      ),
      new Vector3D(this.x_func(s, t), this.y_func(s, t), this.z_func(s, t))
    );

    let ddv_right = Vector3D.difference(dv_right, dv_left);
    let ddv_up = Vector3D.difference(dv_up, dv_down);

    let curving = Math.max(ddv_right.length(), ddv_up.length());

    return curving;
  }
};

export { Vector3D, Vector2D, Matrix, Parametric2D };
