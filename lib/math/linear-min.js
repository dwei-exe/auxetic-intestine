const subtensor_builder = function (current_address, dimensions, builder_function) {
  if (current_address.length === dimensions.length || dimensions[current_address.length] === 0) {
    return builder_function(current_address);
  } else {
    let subtensor = [];
    for (let i_r = 0; i_r < dimensions[current_address.length]; i_r++) {
      subtensor[i_r] = subtensor_builder(current_address.concat([i_r]), dimensions, builder_function);
    }
    return subtensor;
  }
};
const traverse_subtensor = function (current_subtensor, current_address, rank, handler) {
  if (current_address.length === rank) {
    handler(current_subtensor, current_address);
  } else {
    for (let i_r = 0; i_r < current_subtensor.length; i_r++) {
      traverse_subtensor(current_subtensor[i_r], current_address.concat([i_r]), rank, handler);
    }
  }
};

const TENSOR = {};
TENSOR.tensor_builder = function (dimensions, builder_function = (address) => 0) {
  return subtensor_builder([], dimensions, builder_function);
};
TENSOR.traverse_tensor = function (tensor, handler) {
  let rank = TENSOR.rank(tensor);
  traverse_subtensor(tensor, [], rank, handler);
};
TENSOR.coordinate_from_address = function (tensor, address) {
  let current_tensor = tensor;
  for (let i = 0; i < address.length; i++) {
    current_tensor = current_tensor[address[i]];
  }
  return current_tensor;
};
TENSOR.internal_sum = function (tensor) {
  let sum = 0;
  TENSOR.traverse_tensor(tensor, function (coordinate, address) {
    sum += coordinate;
  });
  return sum;
};
TENSOR.rank = function (tensor) {
  let rank = 0;
  let current_tensor = tensor;
  while (Array.isArray(current_tensor)) {
    rank++;
    current_tensor = current_tensor[0];
  }
  return rank;
};
TENSOR.dimensions = function (tensor) {
  let dimensions = [];
  let current_tensor = tensor;
  while (Array.isArray(current_tensor)) {
    dimensions.push(current_tensor.length);
    current_tensor = current_tensor[0];
  }
  return dimensions;
};
TENSOR.add = function (tensor_a, tensor_b) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor_a), function (address) {
    return TENSOR.coordinate_from_address(tensor_a, address) + TENSOR.coordinate_from_address(tensor_b, address);
  });
};
TENSOR.subtract = function (tensor_a, tensor_b) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor_a), function (address) {
    return TENSOR.coordinate_from_address(tensor_a, address) - TENSOR.coordinate_from_address(tensor_b, address);
  });
};
TENSOR.hadamard = function (tensor_a, tensor_b) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor_a), function (address) {
    return TENSOR.coordinate_from_address(tensor_a, address) * TENSOR.coordinate_from_address(tensor_b, address);
  });
};
TENSOR.scalar_mult = function (tensor, scalar) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor), function (address) {
    return TENSOR.coordinate_from_address(tensor, address) * scalar;
  });
};
TENSOR.product = function (tensor_a, tensor_b) {
  let dimensions_a = TENSOR.dimensions(tensor_a);
  let dimensions_b = TENSOR.dimensions(tensor_b);
  return TENSOR.tensor_builder(dimensions_a.concat(dimensions_b), function (address) {
    return (
      TENSOR.coordinate_from_address(tensor_a, address.slice(0, dimensions_a.length)) *
      TENSOR.coordinate_from_address(tensor_b, address.slice(dimensions_b.length))
    );
  });
};
TENSOR.contraction = function (tensor, i_a, i_b) {
  let dimensions = TENSOR.dimensions(tensor);
  //assuming i_a < i_b, and dim i_a = dim i_b
  return TENSOR.tensor_builder(
    dimensions
      .slice(0, i_a)
      .concat(dimensions.slice(i_a + 1, i_b))
      .concat(dimensions.slice(i_b + 1)),
    function (address) {
      let address_beginning = address.slice(0, i_a);
      let address_middle = address.slice(i_a, i_b - 1);
      let address_end = address.slice(i_b - 1);
      let sum = 0;
      for (let k = 0; k < dimensions[i_a]; k++) {
        sum += TENSOR.coordinate_from_address(tensor, address_beginning.concat([k]).concat(address_middle).concat([k]).concat(address_end));
      }
      return sum;
    }
  );
};
TENSOR.apply_map = function (tensor, map) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor), function (address) {
    return map(TENSOR.coordinate_from_address(tensor, address), address);
  });
};
TENSOR.copy = function (tensor) {
  return TENSOR.tensor_builder(TENSOR.dimensions(tensor), function (address) {
    return TENSOR.coordinate_from_address(tensor, address);
  });
};
TENSOR.dot = function (tensor_a, tensor_b) {
  let sum = 0;
  TENSOR.traverse_tensor(tensor_a, function (coordinate, address) {
    sum += coordinate * TENSOR.coordinate_from_address(tensor_b, address);
  });
  return sum;
};
TENSOR.subtensor = function (tensor, ranges) {
  //range: [lower, length]
  let new_dimensions = [];
  for (let i = 0; i < ranges.length; i++) {
    new_dimensions[i] = ranges[i][1];
  }
  return TENSOR.tensor_builder(new_dimensions, function (address) {
    let transposed_address = [];
    for (let i = 0; i < address.length; i++) {
      transposed_address[i] = address[i] + ranges[i][0];
    }
    return TENSOR.coordinate_from_address(tensor, transposed_address);
  });
};

const QUICK_MATRIX = {};
QUICK_MATRIX.matrix_builder = function (height, width, builder_function = (i, j) => 0) {
  let matrix = [];
  for (let i = 0; i < height; i++) {
    matrix[i] = [];
    for (let j = 0; j < width; j++) {
      matrix[i][j] = builder_function(i, j);
    }
  }
  return matrix;
};
QUICK_MATRIX.matrix_mult_vector = function (matrix, vector) {
  let output_vector = [];
  for (let i = 0; i < matrix.length; i++) {
    output_vector[i] = 0;
    for (let j = 0; j < matrix[i].length; j++) {
      output_vector[i] += matrix[i][j] * vector[j];
    }
  }
  return output_vector;
};
QUICK_MATRIX.covector_mult_matrix = function (vector, matrix) {
  let output_vector = [];
  for (let j = 0; j < matrix[0].length; j++) {
    output_vector[j] = 0;
    for (let i = 0; i < matrix.length; i++) {
      output_vector[j] += matrix[i][j] * vector[i];
    }
  }
  return output_vector;
};
QUICK_MATRIX.dot = function (matrix_a, matrix_b) {
  let sum = 0;
  for (let i = 0; i < matrix_a.length; i++) {
    for (let j = 0; j < matrix_a[i].length; j++) {
      sum += matrix_a[i][j] * matrix_b[i][j];
    }
  }
  return sum;
};
QUICK_MATRIX.add = function (matrix_a, matrix_b) {
  return QUICK_MATRIX.matrix_builder(matrix_a.length, matrix_a[0].length, function (i, j) {
    return matrix_a[i][j] + matrix_b[i][j];
  });
};
QUICK_MATRIX.subtract = function (matrix_a, matrix_b) {
  return QUICK_MATRIX.matrix_builder(matrix_a.length, matrix_a[0].length, function (i, j) {
    return matrix_a[i][j] - matrix_b[i][j];
  });
};
QUICK_MATRIX.hadamard = function (matrix_a, matrix_b) {
  return QUICK_MATRIX.matrix_builder(matrix_a.length, matrix_a[0].length, function (i, j) {
    return matrix_a[i][j] * matrix_b[i][j];
  });
};
QUICK_MATRIX.scalar_mult = function (matrix, scalar) {
  return QUICK_MATRIX.matrix_builder(matrix.length, matrix[0].length, function (i, j) {
    return matrix[i][j] * scalar;
  });
};
QUICK_MATRIX.rotate180 = function (matrix) {
  return QUICK_MATRIX.matrix_builder(matrix.length, matrix[0].length, function (i, j) {
    return matrix[matrix.length - 1 - i][matrix[0].length - 1 - j];
  });
};
QUICK_MATRIX.padded_matrix = function (matrix, padding = []) {
  //padding array is from outer to inner values
  let padded_height = matrix.length + 2 * padding.length;
  let padded_width = matrix[0].length + 2 * padding.length;
  let padded_matrix = QUICK_MATRIX.matrix_builder(padded_height, padded_width, function (i, j) {
    let distance_to_border = Math.min(i, j, padded_height - 1 - i, padded_width - 1 - j);
    if (distance_to_border < padding.length) {
      return padding[distance_to_border];
    } else {
      return matrix[i - padding.length][j - padding.length];
    }
  });
  return padded_matrix;
};
QUICK_MATRIX.submatrix = function (matrix, range_i_start_end, range_j_start_end) {
  let sub_matrix = [];
  for (let i = 0; i < range_i_start_end[1] - range_i_start_end[0]; i++) {
    sub_matrix[i] = [];
    for (let j = 0; j < range_j_start_end[1] - range_j_start_end[0]; j++) {
      sub_matrix[i][j] = matrix[range_i_start_end[0] + i][range_j_start_end[0] + j];
    }
  }
  return sub_matrix;
};
QUICK_MATRIX.convolve = function (padded_matrix, kernel, stride = { vertical: 1, horizontal: 1 }) {
  //make sure that the stride divides the sizes of the padded matrix
  let convolved_matrix = QUICK_MATRIX.matrix_builder(
    (padded_matrix.length - (kernel.length - 1)) / stride.vertical,
    (padded_matrix[0].length - (kernel[0].length - 1)) / stride.horizontal,
    function (i, j) {
      return QUICK_MATRIX.dot(
        kernel,
        QUICK_MATRIX.submatrix(
          padded_matrix,
          [i * stride.vertical, i * stride.vertical + kernel.length],
          [j * stride.horizontal, j * stride.horizontal + kernel[0].length]
        )
      );
    }
  );
  return convolved_matrix;
};
QUICK_MATRIX.apply_map = function (matrix, map) {
  return QUICK_MATRIX.matrix_builder(matrix.length, matrix[0].length, function (i, j) {
    return map(matrix[i][j], i, j);
  });
};
QUICK_MATRIX.internal_sum = function (matrix) {
  let sum = 0;
  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[i].length; j++) {
      sum += matrix[i][j];
    }
  }
  return sum;
};
QUICK_MATRIX.determinant = function (matrix) {
  if (matrix.length === 2) {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  } else if (matrix.length === 1) {
    return matrix[0][0];
  } else {
    let sum = 0;
    for (let i = 0; i < matrix.length; i++) {
      if (matrix[i][0] !== 0) {
        let submatrix = [];
        for (let a = 0; a < matrix.length; a++) {
          if (a !== i) {
            submatrix.push(matrix[a].slice(1));
          }
        }
        sum += (-1) ** i * matrix[i][0] * QUICK_MATRIX.determinant(submatrix);
      }
    }
    return sum;
  }
};
QUICK_MATRIX.sub_column_for_vector = function (matrix, vector, column) {
  return QUICK_MATRIX.matrix_builder(matrix.length, matrix[0].length, function (i, j) {
    if (j === column) {
      return vector[i];
    } else {
      return matrix[i][j];
    }
  });
};

const QUICK_VECTOR = {};
QUICK_VECTOR.vector_builder = function (dimension, builder_function) {
  let vector = [];
  for (let i = 0; i < dimension; i++) {
    vector[i] = builder_function(i);
  }
  return vector;
};
QUICK_VECTOR.subtract = function (vector_a, vector_b) {
  let difference = [];
  for (let i = 0; i < vector_a.length; i++) {
    difference[i] = vector_a[i] - vector_b[i];
  }
  return difference;
};
QUICK_VECTOR.add = function (vector_a, vector_b) {
  let sum = [];
  for (let i = 0; i < vector_a.length; i++) {
    sum[i] = vector_a[i] + vector_b[i];
  }
  return sum;
};
QUICK_VECTOR.scalar_mult = function (vector, scalar) {
  let mult = [];
  for (let i = 0; i < vector.length; i++) {
    mult[i] = vector[i] * scalar;
  }
  return mult;
};
QUICK_VECTOR.apply_map = function (vector, map) {
  return QUICK_VECTOR.vector_builder(vector.length, function (i) {
    return map(vector[i], i);
  });
};
QUICK_VECTOR.hadamard = function (vector_a, vector_b) {
  let product = [];
  for (let i = 0; i < vector_a.length; i++) {
    product[i] = vector_a[i] * vector_b[i];
  }
  return product;
};
QUICK_VECTOR.distance = function (vector) {
  let sum = 0;
  for (let i = 0; i < vector.length; i++) {
    sum += vector[i] ** 2;
  }
  return Math.sqrt(sum);
};
QUICK_VECTOR.dot = function (vector_a, vector_b) {
  let sum = 0;
  for (let i = 0; i < vector_a.length; i++) {
    sum += vector_a[i] * vector_b[i];
  }
  return sum;
};
QUICK_VECTOR.projection = function (vector_a, vector_b) {
  return QUICK_VECTOR.scalar_mult(vector_b, QUICK_VECTOR.dot(vector_a, vector_b) / QUICK_VECTOR.distance(vector_b) ** 2);
};

export { TENSOR, QUICK_MATRIX, QUICK_VECTOR };
