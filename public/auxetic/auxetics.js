import * as THREE from "three";
import Eng from "/public/auxetic/eng.js";
import { QUICK_MATRIX, QUICK_VECTOR } from "/lib/math/linear-min.js";

export const mesh_to_graph = function (mesh_file) {
  const vertices = [];

  const edges = [];

  while (true) {
    const line = Eng.Readline(mesh_file);

    if (line === false) break;

    if (line[0] === "#") continue;

    if (line[0] === "v") {
      const coords = [];

      let localNum = "";

      for (let char of line.slice(2)) {
        if (char === " " || char == "\n") {
          coords.push(parseFloat(localNum));
          localNum = "";
        } else {
          localNum += char;
        }
      }

      const vertex = [coords[0], coords[1], coords[2]];

      vertices.push(vertex);
    }

    if (line[0] === "s") continue;

    // Some edges are counted twice
    if (line[0] === "f") {
      let indices = [];

      let localNum = "";

      for (let char of line.slice(2)) {
        if (char === " " || char == "\n") {
          indices.push(parseInt(localNum) - 1);
          localNum = "";
        } else {
          localNum += char;
        }
      }

      for (let i = 0; i < indices.length - 1; i++) {
        edges.push([indices[i], indices[i + 1]]);
      }
      edges.push([indices[0], indices[indices.length - 1]]);
    }
  }

  return { vertices: vertices, edges: edges };
};

const non_repeating_edges = function (edges) {
  const surviving_indices = new Array(edges.length).fill(true);
  for (let i = 0; i < edges.length; i++) {
    for (let j = 0; j < edges.length; j++) {
      if (i !== j) {
        if (edges[i].includes(edges[j][0]) && edges[i].includes(edges[j][1])) {
          surviving_indices[i] = false;
          surviving_indices[j] = false;
        }
      }
    }
  }
  return surviving_indices;
};

const edge_length_3D = function (edge, vertices) {
  const vertex1 = vertices[edge[0]];
  const vertex2 = vertices[edge[1]];
  return Math.sqrt((vertex1[0] - vertex2[0]) ** 2 + (vertex1[1] - vertex2[1]) ** 2 + (vertex1[2] - vertex2[2]) ** 2);
};

const neighbors = function (vertex_index, graph) {
  const edges = graph.edges;

  const neighbors_indices = [];
  for (const edge of edges) {
    if (edge.includes(vertex_index) === true) {
      const other_vertex_index = edge[0] === vertex_index ? edge[1] : edge[0];
      if (neighbors_indices.includes(other_vertex_index) === false) {
        neighbors_indices.push(other_vertex_index);
      }
    }
  }

  return neighbors_indices;
};

const edge_is_in_array = function (array, edge) {
  for (const a of array) {
    if (a.includes(edge[0]) && a.includes(edge[1])) {
      return true;
    }
  }
  return false;
};

export const boundary = function (graph, width, height, starting_vertex_index) {
  const vertices = graph.vertices; // actual coordinates in 3D [v1, v2, ..., vn]
  const edges = graph.edges; // [ [index of vi, index of vj] ]

  const boundary_edges = [];
  const boundary_edges_boolean = non_repeating_edges(edges);
  for (let i = 0; i < boundary_edges_boolean.length; i++) {
    if (boundary_edges_boolean[i] === true) {
      boundary_edges.push([edges[i][0], edges[i][1]]);
    }
  }

  const boundary_vertices_indices = []; // indices
  for (const boundary_edge of boundary_edges) {
    if (boundary_vertices_indices.includes(boundary_edge[0]) === false) {
      boundary_vertices_indices.push(boundary_edge[0]);
    }
    if (boundary_vertices_indices.includes(boundary_edge[1]) === false) {
      boundary_vertices_indices.push(boundary_edge[1]);
    }
  }

  const ordered_boundary_vertices_indices = []; // indices still with respect to "vertices", but these indices are in new order
  ordered_boundary_vertices_indices.push(boundary_vertices_indices[starting_vertex_index]);
  while (ordered_boundary_vertices_indices.length < boundary_vertices_indices.length) {
    const last_vertex_index = ordered_boundary_vertices_indices[ordered_boundary_vertices_indices.length - 1];

    for (const neighbor_vertex_index of neighbors(last_vertex_index, graph)) {
      if (
        edge_is_in_array(boundary_edges, [last_vertex_index, neighbor_vertex_index]) &&
        ordered_boundary_vertices_indices.includes(neighbor_vertex_index) === false
      ) {
        ordered_boundary_vertices_indices.push(neighbor_vertex_index);
        break;
      }
    }
  }

  const boundary_edge_lengths = [];
  let boundary_length = 0;
  for (let i = 0; i < ordered_boundary_vertices_indices.length; i++) {
    const vertex1_index = ordered_boundary_vertices_indices[i];
    const vertex2_index =
      i + 1 < ordered_boundary_vertices_indices.length ? ordered_boundary_vertices_indices[i + 1] : ordered_boundary_vertices_indices[0];
    const edge_length = edge_length_3D([vertex1_index, vertex2_index], vertices);
    boundary_edge_lengths.push(edge_length);
    boundary_length += edge_length;
  }

  const parametrized_ordered_boundary_vertices = []; // 1-to-1 with ordered_boundary_vertices_indices
  const turtle = [0, 0]; // positions X and Y
  let counter = 0;
  for (let i = 0; i < width - 1; i++) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[0] += boundary_edge_lengths[counter];
    counter++;
  }
  for (let i = 0; i < height - 1; i++) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[1] += boundary_edge_lengths[counter];
    counter++;
  }
  for (let i = 0; i < width - 1; i++) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[0] -= boundary_edge_lengths[counter];
    counter++;
  }
  for (let i = 0; i < height; i++) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[1] -= boundary_edge_lengths[counter];
    counter++;
  }

  return {
    ordered_boundary_vertices_indices: ordered_boundary_vertices_indices, // points to which vertex in the original array the parametrized vertex is related to
    parametrized_ordered_boundary_vertices: parametrized_ordered_boundary_vertices,
  };
};

export const parametrization = function (graph, boundary, metric) {
  const vertices = graph.vertices;
  const edges = graph.edges;

  const ordered_boundary_vertices_indices = boundary.ordered_boundary_vertices_indices;
  const parametrized_ordered_boundary_vertices = boundary.parametrized_ordered_boundary_vertices;

  const interior_vertices_indices = []; // indices
  for (let i = 0; i < vertices.length; i++) {
    if (ordered_boundary_vertices_indices.includes(i) === false) {
      interior_vertices_indices.push(i);
    }
  }

  const sorted_edges = QUICK_MATRIX.matrix_builder(vertices.length, vertices.length, function (i, j) {
    for (const edge of edges) {
      if (edge.includes(i) && edge.includes(j)) {
        return edge;
      }
    }
    return false;
  }); // with respect to vertices, actually a rank-3 tensor, matrix of edges

  // change spring constant D as needed
  const spring_constants = QUICK_MATRIX.apply_map(sorted_edges, function (edge, i, j) {
    if (edge !== false) {
      const edge_length = edge_length_3D(edge, vertices);
      const D = metric(edge_length);
      return D;
    } else {
      return false;
    }
  }); // with respect to vertices

  const weights = QUICK_MATRIX.matrix_builder(vertices.length, vertices.length, function (i, j) {
    let total_neighboring_spring_constants = 0;
    for (const neighbor_vertex_index of neighbors(i, graph)) {
      total_neighboring_spring_constants += spring_constants[i][neighbor_vertex_index];
    }
    const lambda = spring_constants[i][j] / total_neighboring_spring_constants;
    return lambda;
  }); // with respect to vertices

  const filtered_weights = QUICK_MATRIX.matrix_builder(interior_vertices_indices.length, interior_vertices_indices.length, function (i, j) {
    if (i === j) {
      return 1;
    } else if (neighbors(interior_vertices_indices[i], graph).includes(interior_vertices_indices[j])) {
      return -weights[interior_vertices_indices[i]][interior_vertices_indices[j]];
    } else {
      return 0;
    }
  }); // with respect to interior vertices

  const boundary_total_weights_u = QUICK_VECTOR.apply_map(interior_vertices_indices, function (interior_vertex_index, i) {
    let accumulation = 0;
    for (const neighbor_vertex_index of neighbors(interior_vertex_index, graph)) {
      if (ordered_boundary_vertices_indices.includes(neighbor_vertex_index) === true) {
        accumulation +=
          weights[interior_vertex_index][neighbor_vertex_index] *
          parametrized_ordered_boundary_vertices[ordered_boundary_vertices_indices.indexOf(neighbor_vertex_index)][0];
      }
    }
    return accumulation;
  }); // same order as interior vertices

  const boundary_total_weights_v = QUICK_VECTOR.apply_map(interior_vertices_indices, function (interior_vertex_index, i) {
    let accumulation = 0;
    for (const neighbor_vertex_index of neighbors(interior_vertex_index, graph)) {
      if (ordered_boundary_vertices_indices.includes(neighbor_vertex_index) === true) {
        accumulation +=
          weights[interior_vertex_index][neighbor_vertex_index] *
          parametrized_ordered_boundary_vertices[ordered_boundary_vertices_indices.indexOf(neighbor_vertex_index)][1];
      }
    }
    return accumulation;
  }); // same order as interior vertices

  const parametrized_interior_vertices_u = math.lusolve(filtered_weights, boundary_total_weights_u);
  const parametrized_interior_vertices_v = math.lusolve(filtered_weights, boundary_total_weights_v);

  const parametrized_interior_vertices = []; // same order as interior vertices
  for (let i = 0; i < parametrized_interior_vertices_u.length; i++) {
    parametrized_interior_vertices.push([parametrized_interior_vertices_u[i][0], parametrized_interior_vertices_v[i][0]]);
  }

  const indices = interior_vertices_indices.concat(ordered_boundary_vertices_indices);
  const parametrized_vertices = parametrized_interior_vertices.concat(parametrized_ordered_boundary_vertices);

  return {
    vertices: vertices, // original vertices
    edges: edges, // edges don't change
    indices: indices, // the indices tell us the original vertex that the parametrized vertex is mapped to
    parametrized_vertices: parametrized_vertices, // parametrized vertices
  };
};

const edge_length_2D = function (edge, vertices) {
  const vertex1 = vertices[edge[0]];
  const vertex2 = vertices[edge[1]];
  return Math.sqrt((vertex1[0] - vertex2[0]) ** 2 + (vertex1[1] - vertex2[1]) ** 2);
};

export const auxetic_info = function (parametrized, angle_restriction) {
  const vertices = parametrized.vertices;
  const edges = parametrized.edges;
  const parametrized_vertices_indices = parametrized.indices; // go from parametrized to original
  const parametrized_vertices = parametrized.parametrized_vertices;

  const graph = {
    vertices: vertices,
    edges: edges,
  };

  // go from original to parametrized
  const inverse_indices = new Array(parametrized_vertices_indices.length);
  for (let i = 0; i < parametrized_vertices_indices.length; i++) {
    inverse_indices[parametrized_vertices_indices[i]] = i;
  }

  const parametrized_edges = [];
  for (const edge of edges) {
    const new_vertex1_index = inverse_indices[edge[0]];
    const new_vertex2_index = inverse_indices[edge[1]];
    parametrized_edges.push([new_vertex1_index, new_vertex2_index]);
  }

  const original_edge_lengths = QUICK_VECTOR.apply_map(edges, function (original_edge, i, j) {
    return edge_length_3D(original_edge, vertices);
  });

  const new_edge_lengths = QUICK_VECTOR.apply_map(edges, function (original_edge, i, j) {
    const new_edge = [inverse_indices[original_edge[0]], inverse_indices[original_edge[1]]];
    return edge_length_2D(new_edge, parametrized_vertices);
  });

  // edge lengths are 1-to-1 indices with each other

  const edge_length_difference = QUICK_VECTOR.subtract(original_edge_lengths, new_edge_lengths); // from flat to 3D

  //const max_rotation_angles = Math.PI / 5; // change later
  const max_rotation_angles = QUICK_VECTOR.apply_map(parametrized_vertices, function (parametrized_vertex, i) {
    const parametrized_vertex_neighbor_indices = neighbors(i, { edges: parametrized_edges });
    const angles = new Array(parametrized_vertices.length);
    for (let j = 0; j < parametrized_vertex_neighbor_indices.length; j++) {
      const parametrized_neighbor_vertex = parametrized_vertices[parametrized_vertex_neighbor_indices[j]];
      angles[parametrized_vertex_neighbor_indices[j]] = Math.atan2(
        parametrized_neighbor_vertex[1] - parametrized_vertex[1],
        parametrized_neighbor_vertex[0] - parametrized_vertex[0]
      );
    }
    parametrized_vertex_neighbor_indices.sort(function (a, b) {
      return angles[a] - angles[b];
    });
    const differences = [];
    for (let j = 0; j < parametrized_vertex_neighbor_indices.length; j++) {
      let diff = 0;
      if (j + 1 < parametrized_vertex_neighbor_indices.length) {
        diff = angles[parametrized_vertex_neighbor_indices[j + 1]] - angles[parametrized_vertex_neighbor_indices[j]];
      } else {
        diff = angles[parametrized_vertex_neighbor_indices[0]] - angles[parametrized_vertex_neighbor_indices[j]] + Math.PI * 2;
      }
      diff *= angle_restriction;
      differences.push(diff);
    }

    return Math.min(...differences, Math.PI / 3);
  });

  const max_radii = QUICK_VECTOR.apply_map(parametrized_vertices, function (parametrized_vertex, i) {
    const parametrized_vertex_neighbor_indices = neighbors(i, { edges: parametrized_edges });
    const lengths = [];
    for (const parametrized_vertex_neighbor_index of parametrized_vertex_neighbor_indices) {
      lengths.push(edge_length_2D([i, parametrized_vertex_neighbor_index], parametrized_vertices));
    }
    const min_length = Math.min(...lengths);
    return min_length / 3;
  });

  const parametrized_edge_length_differences = [];
  const radii = QUICK_VECTOR.apply_map(parametrized_edges, function (parametrized_edge, i) {
    parametrized_edge_length_differences.push(
      edge_length_3D([parametrized_vertices_indices[parametrized_edge[0]], parametrized_vertices_indices[parametrized_edge[1]]], vertices) -
        edge_length_2D(parametrized_edge, parametrized_vertices)
    );

    const E2 = edge_length_2D(parametrized_edge, parametrized_vertices) / 2;

    const r1 =
      (E2 ** 2 - (edge_length_difference[i] + E2) ** 2) /
      (2 * E2 * Math.cos(max_rotation_angles[parametrized_edge[0]]) - 2 * (edge_length_difference[i] + E2));
    const r2 =
      (E2 ** 2 - (edge_length_difference[i] + E2) ** 2) /
      (2 * E2 * Math.cos(max_rotation_angles[parametrized_edge[1]]) - 2 * (edge_length_difference[i] + E2));
    return [Math.max(Math.min(r1, max_radii[parametrized_edge[0]]), 0), Math.max(Math.min(r2, max_radii[parametrized_edge[1]]), 0)];
  });

  return {
    parametrized_vertices: parametrized_vertices,
    parametrized_edges: parametrized_edges,
    max_radii: max_radii,
    radii: radii,
    max_rotation_angles: max_rotation_angles,
    parametrized_edge_length_differences: parametrized_edge_length_differences,
  };
};

export const auxetic_mesh = function (scene, models, auxetic_infoed, thickness) {
  const parametrized_vertices = auxetic_infoed.parametrized_vertices;
  const parametrized_edges = auxetic_infoed.parametrized_edges;
  const radii = auxetic_infoed.radii;
  const max_radii = auxetic_infoed.max_radii;
  const max_rotation_angles = auxetic_infoed.max_rotation_angles;

  const unit = models.unit;
  const knob = models.knob;
  const rod_ext = models.rod_ext;
  const rod_mid = models.rod_mid;

  for (let i = 0; i < parametrized_vertices.length; i++) {
    const parametrized_vertex = parametrized_vertices[i];

    const cloned_unit = unit.clone();
    cloned_unit.translateX(parametrized_vertex[0]);
    cloned_unit.translateY(parametrized_vertex[1]);
    cloned_unit.scale.set(thickness / 50, thickness / 50, thickness / 50);
    scene.add(cloned_unit);

    const neighbor_vertex_indices = neighbors(i, { edges: parametrized_edges });
    for (const neighbor_vertex_index of neighbor_vertex_indices) {
      const neighbor_vertex = parametrized_vertices[neighbor_vertex_index];
      const delta_x = neighbor_vertex[0] - parametrized_vertex[0];
      const delta_y = neighbor_vertex[1] - parametrized_vertex[1];
      let angle = Math.atan2(delta_y, delta_x);

      let wanted_edge_index;
      let wanted_edgy_index;
      for (let j = 0; j < parametrized_edges.length; j++) {
        if (parametrized_edges[j].includes(i) === true && parametrized_edges[j].includes(neighbor_vertex_index) === true) {
          wanted_edgy_index = parametrized_edges[j][0] === i ? 0 : 1;
          wanted_edge_index = j;
          break;
        }
      }

      const radius = radii[wanted_edge_index][wanted_edgy_index];

      if (radius >= 0) {
        angle += max_rotation_angles[i];
      }

      const unit_rod_wrapper = new THREE.Object3D();
      unit_rod_wrapper.translateX(parametrized_vertex[0]);
      unit_rod_wrapper.translateY(parametrized_vertex[1]);
      unit_rod_wrapper.rotation.set(0, 0, angle - Math.PI / 2, "YXZ");

      const cloned_rod_node = rod_mid.clone();
      cloned_rod_node.scale.set(Math.max(radius * 0.75, thickness / 150), radius * 1.1, thickness / 10);

      const cloned_knob = knob.clone();
      cloned_knob.translateY(radius * 1.25);
      cloned_knob.scale.set(Math.max(radius * 0.75, thickness / 150), radius * 0.9, thickness / 60);

      unit_rod_wrapper.add(cloned_rod_node);
      unit_rod_wrapper.add(cloned_knob);
      scene.add(unit_rod_wrapper);

      const average_x = (parametrized_vertex[0] + neighbor_vertex[0]) / 2;
      const average_y = (parametrized_vertex[1] + neighbor_vertex[1]) / 2;

      const base_x = parametrized_vertex[0] + 1.24 * radius * Math.cos(angle);
      const base_y = parametrized_vertex[1] + 1.24 * radius * Math.sin(angle);

      const rod_length = Math.sqrt((average_x - base_x) ** 2 + (average_y - base_y) ** 2);

      const rod_ext_wrapper = new THREE.Object3D();
      rod_ext_wrapper.translateX(base_x);
      rod_ext_wrapper.translateY(base_y);
      rod_ext_wrapper.rotation.set(0, 0, Math.atan2(average_y - base_y, average_x - base_x) - Math.PI / 2);

      const cloned_rod_ext = rod_ext.clone();
      cloned_rod_ext.scale.set(Math.max(radius * 0.7, thickness / 150), radius * 0.65, thickness / 30);
      cloned_rod_ext.translateY(0.005 * radius);

      const cloned_rod_mid = rod_mid.clone();
      cloned_rod_mid.translateY(0.15 * radius);
      cloned_rod_mid.scale.set(Math.max(radius * 0.5, thickness / 150), rod_length * 1.1, thickness / 30);

      rod_ext_wrapper.add(cloned_rod_ext);
      rod_ext_wrapper.add(cloned_rod_mid);
      scene.add(rod_ext_wrapper);
    }
  }
};

export const folded_mesh = function (scene, models, parametrized, auxetic_infoed, thickness) {
  const vertices = parametrized.vertices;
  const indices = parametrized.indices;
  const parametrized_vertices = auxetic_infoed.parametrized_vertices;
  const parametrized_edges = auxetic_infoed.parametrized_edges;
  const radii = auxetic_infoed.radii;

  const unit = models.unit;
  const knob = models.knob;
  const rod_ext = models.rod_ext;
  const rod_mid = models.rod_mid;

  const return_normal = [];

  const folded_models = [];

  for (let i = 0; i < parametrized_vertices.length; i++) {
    const folded_model = {};

    const vertex = vertices[indices[i]];

    const cloned_unit = unit.clone();
    scene.add(cloned_unit);

    folded_model.unit = cloned_unit;

    const neighbor_vertex_indices = neighbors(i, { edges: parametrized_edges });

    const normals = [];
    let current_vector;

    for (const neighbor_vertex_index of neighbor_vertex_indices) {
      const neighbor_vertex = vertices[indices[neighbor_vertex_index]];

      const delta_x = neighbor_vertex[0] - vertex[0];
      const delta_y = neighbor_vertex[1] - vertex[1];
      const delta_z = neighbor_vertex[2] - vertex[2];

      if (current_vector !== undefined) {
        normals.push([
          delta_y * current_vector[2] - delta_z * current_vector[1],
          delta_z * current_vector[0] - delta_x * current_vector[2],
          delta_x * current_vector[1] - delta_y * current_vector[0],
        ]);
      }
      current_vector = [delta_x, delta_y, delta_z];
    }

    let average_normal = [0, 0, 0];
    for (const normal of normals) {
      average_normal = QUICK_VECTOR.add(average_normal, normal);
    }
    average_normal = QUICK_VECTOR.scalar_mult(average_normal, 1 / normals.length);

    return_normal.push(average_normal);

    const theta1 = Math.atan2(average_normal[2], average_normal[0]);
    const theta2 = Math.atan2(average_normal[1], Math.sqrt(average_normal[0] ** 2 + average_normal[2] ** 2));

    cloned_unit.translateX(vertex[0]);
    cloned_unit.translateY(vertex[1]);
    cloned_unit.translateZ(vertex[2]);
    cloned_unit.rotation.set(-theta2, Math.PI / 2 - theta1, 0, "YXZ");
    cloned_unit.scale.set(thickness / 50, thickness / 50, thickness / 50);

    const normalized_normal = QUICK_VECTOR.scalar_mult(average_normal, 1 / QUICK_VECTOR.distance(average_normal));

    folded_model.normal_axis = normalized_normal;

    folded_model.unit_rods = [];
    folded_model.rod_exts = [];
    for (const neighbor_vertex_index of neighbor_vertex_indices) {
      const neighbor_vertex = vertices[indices[neighbor_vertex_index]];

      const delta_x = neighbor_vertex[0] - vertex[0];
      const delta_y = neighbor_vertex[1] - vertex[1];
      const delta_z = neighbor_vertex[2] - vertex[2];
      let angle1 = Math.atan2(delta_z, delta_x);
      let angle2 = Math.atan2(delta_y, Math.sqrt(delta_x ** 2 + delta_z ** 2));

      let wanted_edge_index;
      let wanted_edgy_index;
      for (let j = 0; j < parametrized_edges.length; j++) {
        if (parametrized_edges[j].includes(i) === true && parametrized_edges[j].includes(neighbor_vertex_index) === true) {
          wanted_edgy_index = parametrized_edges[j][0] === i ? 0 : 1;
          wanted_edge_index = j;
          break;
        }
      }

      const radius = radii[wanted_edge_index][wanted_edgy_index];

      if (radius >= 0) {
        //angle1 += max_rotation_angles[i];
      }

      const unit_rod_wrapper = new THREE.Object3D();
      unit_rod_wrapper.translateX(vertex[0]);
      unit_rod_wrapper.translateY(vertex[1]);
      unit_rod_wrapper.translateZ(vertex[2]);
      unit_rod_wrapper.lookAt(
        new THREE.Vector3(normalized_normal[0] + vertex[0], normalized_normal[2] + vertex[2], normalized_normal[2] + vertex[2])
      );
      //unit_rod_wrapper.rotation.set(Math.PI / 2 - angle2, Math.PI / 2 - angle1, 0, "YXZ");

      const cloned_rod_node = rod_mid.clone();
      cloned_rod_node.scale.set(Math.max(radius * 0.75, thickness / 150), radius * 1.1, thickness / 10);

      const cloned_knob = knob.clone();
      cloned_knob.translateY(radius * 1.25);
      cloned_knob.scale.set(Math.max(radius * 0.75, thickness / 150), radius * 0.9, thickness / 60);

      unit_rod_wrapper.add(cloned_rod_node);
      unit_rod_wrapper.add(cloned_knob);
      scene.add(unit_rod_wrapper);

      folded_model.unit_rods.push(unit_rod_wrapper);

      const average_x = (vertex[0] + neighbor_vertex[0]) / 2;
      const average_y = (vertex[1] + neighbor_vertex[1]) / 2;
      const average_z = (vertex[2] + neighbor_vertex[2]) / 2;

      const base_x = vertex[0]; // + 1.24 * radius * Math.cos(angle1);
      const base_y = vertex[1]; // + 1.24 * radius * Math.sin(angle2);
      const base_z = vertex[2]; // + 1.24 * radius * Math.sin(angle1);

      const rod_length = Math.sqrt((average_x - base_x) ** 2 + (average_y - base_y) ** 2 + (average_z - base_z) ** 2);

      const phi1 = Math.atan2(average_z - base_z, average_x - base_x);
      const phi2 = Math.atan2(average_y - base_y, Math.sqrt((average_z - base_z) ** 2 + (average_x - base_x) ** 2));

      const rod_ext_wrapper = new THREE.Object3D();
      rod_ext_wrapper.translateX(base_x);
      rod_ext_wrapper.translateY(base_y);
      rod_ext_wrapper.translateZ(base_z);
      rod_ext_wrapper.scale.set(Math.max(radius * 0.5, thickness / 150), rod_length * 1.1, thickness / 30);
      rod_ext_wrapper.rotation.set(Math.PI / 2 - phi2, Math.PI / 2 - phi1, 0, "YXZ");

      const cloned_rod_ext = rod_ext.clone();
      cloned_rod_ext.scale.set(Math.max(radius * 0.7, thickness / 150), radius * 0.65, thickness / 30);
      cloned_rod_ext.translateY(0.005 * radius);

      const cloned_rod_mid = rod_mid.clone();
      cloned_rod_mid.translateY(1.24 * radius);

      rod_ext_wrapper.add(cloned_rod_ext);
      rod_ext_wrapper.add(cloned_rod_mid);
      scene.add(rod_ext_wrapper);

      folded_model.rod_exts.push(rod_ext_wrapper);
    }

    folded_models.push(folded_model);
  }

  return folded_models;
};

export const morph_intestine = function (parametrized, auxetic_infoed, folded_models, time, delta_time, thickness) {
  const vertices = parametrized.vertices;
  const indices = parametrized.indices;
  const parametrized_edges = auxetic_infoed.parametrized_edges;
  const radii = auxetic_infoed.radii;

  const transformed_vertices = [];

  for (let i = 0; i < folded_models.length; i++) {
    const folded_model = folded_models[i];

    const vertex = vertices[indices[i]];

    const cloned_unit = folded_model.unit;

    cloned_unit.translateY(Math.cos(time * 4 + vertex[0] * 20) * delta_time * 0.05);
    cloned_unit.translateZ(Math.cos(time * 4 + vertex[0] * 20) * delta_time * 0.05);

    transformed_vertices.push([cloned_unit.position.x, cloned_unit.position.y, cloned_unit.position.z]);
  }

  //folded_models is indexed with respect to parametrized vertices
  for (let i = 0; i < folded_models.length; i++) {
    const folded_model = folded_models[i];

    const transformed_vertex = transformed_vertices[i];

    const neighbor_vertex_indices = neighbors(i, { edges: parametrized_edges });

    for (let j = 0; j < neighbor_vertex_indices.length; j++) {
      const neighbor_vertex_index = neighbor_vertex_indices[j];
      const neighbor_transformed_vertex = transformed_vertices[neighbor_vertex_index];

      const delta_x = neighbor_transformed_vertex[0] - transformed_vertex[0];
      const delta_y = neighbor_transformed_vertex[1] - transformed_vertex[1];
      const delta_z = neighbor_transformed_vertex[2] - transformed_vertex[2];
      let angle1 = Math.atan2(delta_z, delta_x);
      let angle2 = Math.atan2(delta_y, Math.sqrt(delta_x ** 2 + delta_z ** 2));

      let wanted_edge_index;
      let wanted_edgy_index;
      for (let j = 0; j < parametrized_edges.length; j++) {
        if (parametrized_edges[j].includes(i) === true && parametrized_edges[j].includes(neighbor_vertex_index) === true) {
          wanted_edgy_index = parametrized_edges[j][0] === i ? 0 : 1;
          wanted_edge_index = j;
          break;
        }
      }

      const radius = radii[wanted_edge_index][wanted_edgy_index];

      if (radius >= 0) {
        //angle1 += max_rotation_angles[i];
      }

      const unit_rod_wrapper = folded_model.unit_rods[j];
      unit_rod_wrapper.position.set(transformed_vertex[0], transformed_vertex[1], transformed_vertex[2]);
      unit_rod_wrapper.rotation.set(Math.PI / 2 - angle2, Math.PI / 2 - angle1, 0, "YXZ");
      //unit_rod_wrapper.rotateOnAxis(new THREE.Vector3(folded_model.normal_axis[0], folded_model.normal_axis[1], folded_model.normal_axis[2]), delta_time);

      const average_x = (transformed_vertex[0] + neighbor_transformed_vertex[0]) / 2;
      const average_y = (transformed_vertex[1] + neighbor_transformed_vertex[1]) / 2;
      const average_z = (transformed_vertex[2] + neighbor_transformed_vertex[2]) / 2;

      const base_x = transformed_vertex[0]; // + 1.24 * radius * Math.cos(angle1);
      const base_y = transformed_vertex[1]; // + 1.24 * radius * Math.sin(angle2);
      const base_z = transformed_vertex[2]; // + 1.24 * radius * Math.sin(angle1);

      const rod_length = Math.sqrt((average_x - base_x) ** 2 + (average_y - base_y) ** 2 + (average_z - base_z) ** 2);

      const phi1 = Math.atan2(average_z - base_z, average_x - base_x);
      const phi2 = Math.atan2(average_y - base_y, Math.sqrt((average_z - base_z) ** 2 + (average_x - base_x) ** 2));

      const rod_ext_wrapper = folded_model.rod_exts[j];
      rod_ext_wrapper.position.set(base_x, base_y, base_z);
      rod_ext_wrapper.rotation.set(Math.PI / 2 - phi2, Math.PI / 2 - phi1, 0, "YXZ");
      rod_ext_wrapper.scale.set(Math.max(radius * 0.5, thickness / 150), rod_length * 1.1, thickness / 30);
    }
  }
};

const edge_direction = function (vertices, edge) {
  const first_vertex = vertices[edge[0]];
  const last_vertex = vertices[edge[1]];

  let direction = QUICK_VECTOR.subtract(last_vertex, first_vertex);
  direction = QUICK_VECTOR.scalar_mult(direction, 1 / QUICK_VECTOR.distance(direction));
  return direction;
};

//external force constant (assuming for small movements)
export const morph_forces = function (morphed_graph, parametrized, auxetic_infoed, external_forces, unit_rotations) {
  const vertices = morphed_graph.vertices;
  const edges = morphed_graph.edges;
  const radii = auxetic_infoed.radii;
  const indices = parametrized.indices;

  // go from original to parametrized
  const inverse_indices = new Array(indices.length);
  for (let i = 0; i < indices.length; i++) {
    inverse_indices[indices[i]] = i;
  }

  const weight_matrix = QUICK_MATRIX.matrix_builder(vertices.length * 3, edges.length, function (j, i) {
    const current_edge = edges[i];
    const direction = edge_direction(vertices, current_edge);

    if (current_edge.includes(j)) {
      return direction[j % 3];
    } else {
      return 0;
    }
  });
  const counter_forces = QUICK_VECTOR.vector_builder(vertices.length * 3, function (i) {
    return -external_forces[Math.floor(i / 3)][i % 3];
  });

  const solved_rod_forces_components = math.lusolve(weight_matrix, counter_forces);

  const solved_rod_forces_magnitudes = QUICK_VECTOR.vector_builder(edges.length, function (i) {
    return QUICK_VECTOR.distance([
      solved_rod_forces_components[i * 3],
      solved_rod_forces_components[i * 3 + 1],
      solved_rod_forces_components[i * 3 + 2],
    ]);
  });

  // very approximate torque
  const torques = QUICK_VECTOR.vector_builder(vertices.length, function (i) {
    let sum = 0;
    for (let i = 0; i < edges.length; i++) {
      const edge = edges[i];
      if (edge.includes(i)) {
        const edgy_index = edge[0] === i ? edge[1] : edge[0];
        const radius = radii[inverse_indices[i]][inverse_indices[edgy_index]];
        sum += Math.cos(unit_rotations[i]) * solved_rod_forces_magnitudes[i] * radius;
      }
    }
    return sum;
  });

  return {
    rod_forces_magnitudes: solved_rod_forces_magnitudes,
    torques: torques,
  };
};
