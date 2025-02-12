import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js";
import { OBJLoader } from "three/addons/loaders/OBJLoader.js";
import { GLTFExporter } from "three/addons/exporters/GLTFExporter.js";
import { Sky } from "three/addons/objects/Sky.js";

import Eng from "/public/eng.js";
import { QUICK_MATRIX, QUICK_VECTOR } from "/public/linear-min.js";

let boundary_lengthh;

const mesh_to_graph = function (mesh_file) {
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
          coords.push(localNum);
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
  return Math.sqrt(
    (vertex1[0] - vertex2[0]) ** 2 + (vertex1[1] - vertex2[1]) ** 2 + (vertex1[2] - vertex2[2]) ** 2
  );
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

const boundary = function (graph, topology = "disk", rect_ratio = 1) {
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
  ordered_boundary_vertices_indices.push(boundary_vertices_indices[1]);
  while (ordered_boundary_vertices_indices.length < boundary_vertices_indices.length) {
    const last_vertex_index =
      ordered_boundary_vertices_indices[ordered_boundary_vertices_indices.length - 1];

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
      i + 1 < ordered_boundary_vertices_indices.length
        ? ordered_boundary_vertices_indices[i + 1]
        : ordered_boundary_vertices_indices[0];
    const edge_length = edge_length_3D([vertex1_index, vertex2_index], vertices);
    boundary_edge_lengths.push(edge_length);
    boundary_length += edge_length;
  }
  boundary_lengthh = boundary_length;

  const parametrized_ordered_boundary_vertices = []; // 1-to-1 with ordered_boundary_vertices_indices
  const turtle = [0, 0]; // positions X and Y
  let counter = 0;
  let x_counter = 0;
  let y_counter = 0;
  while (x_counter < 7) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[0] += boundary_edge_lengths[counter];
    x_counter++;
    counter++;
  }
  while (y_counter < 9) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[1] += boundary_edge_lengths[counter];
    y_counter++;
    counter++;
  }
  while (x_counter > 0) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[0] -= boundary_edge_lengths[counter];
    x_counter--;
    counter++;
  }
  while (y_counter > 0 && counter < boundary_edge_lengths.length) {
    parametrized_ordered_boundary_vertices.push([turtle[0], turtle[1]]);
    turtle[1] -= boundary_edge_lengths[counter];
    y_counter--;
    counter++;
  }

  return {
    ordered_boundary_vertices_indices: ordered_boundary_vertices_indices, // points to which vertex in the original array the parametrized vertex is related to
    parametrized_ordered_boundary_vertices: parametrized_ordered_boundary_vertices,
  };
};

const parametrization = function (graph, boundary) {
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

  const sorted_edges = QUICK_MATRIX.matrix_builder(
    vertices.length,
    vertices.length,
    function (i, j) {
      for (const edge of edges) {
        if (edge.includes(i) && edge.includes(j)) {
          return edge;
        }
      }
      return false;
    }
  ); // with respect to vertices, actually a rank-3 tensor, matrix of edges

  // change spring constant D as needed
  const spring_constants = QUICK_MATRIX.apply_map(sorted_edges, function (edge, i, j) {
    if (edge !== false) {
      const edge_length = edge_length_3D(edge, vertices);
      const D = 1 / edge_length;
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

  const filtered_weights = QUICK_MATRIX.matrix_builder(
    interior_vertices_indices.length,
    interior_vertices_indices.length,
    function (i, j) {
      if (i === j) {
        return 1;
      } else if (
        neighbors(interior_vertices_indices[i], graph).includes(interior_vertices_indices[j])
      ) {
        return -weights[interior_vertices_indices[i]][interior_vertices_indices[j]];
      } else {
        return 0;
      }
    }
  ); // with respect to interior vertices

  const boundary_total_weights_u = QUICK_VECTOR.apply_map(
    interior_vertices_indices,
    function (interior_vertex_index, i) {
      let accumulation = 0;
      for (const neighbor_vertex_index of neighbors(interior_vertex_index, graph)) {
        if (ordered_boundary_vertices_indices.includes(neighbor_vertex_index) === true) {
          accumulation +=
            weights[interior_vertex_index][neighbor_vertex_index] *
            parametrized_ordered_boundary_vertices[
              ordered_boundary_vertices_indices.indexOf(neighbor_vertex_index)
            ][0];
        }
      }
      return accumulation;
    }
  ); // same order as interior vertices

  const boundary_total_weights_v = QUICK_VECTOR.apply_map(
    interior_vertices_indices,
    function (interior_vertex_index, i) {
      let accumulation = 0;
      for (const neighbor_vertex_index of neighbors(interior_vertex_index, graph)) {
        if (ordered_boundary_vertices_indices.includes(neighbor_vertex_index) === true) {
          accumulation +=
            weights[interior_vertex_index][neighbor_vertex_index] *
            parametrized_ordered_boundary_vertices[
              ordered_boundary_vertices_indices.indexOf(neighbor_vertex_index)
            ][1];
        }
      }
      return accumulation;
    }
  ); // same order as interior vertices

  const parametrized_interior_vertices_u = math.lusolve(filtered_weights, boundary_total_weights_u);
  const parametrized_interior_vertices_v = math.lusolve(filtered_weights, boundary_total_weights_v);

  /*
  const filtered_weights_determinant = QUICK_MATRIX.determinant(filtered_weights);
  console.table(filtered_weights_determinant);
  const parametrized_interior_vertices_u = QUICK_VECTOR.vector_builder(
    interior_vertices_indices.length,
    function (i) {
      return (
        QUICK_MATRIX.determinant(
          QUICK_MATRIX.sub_column_for_vector(filtered_weights, boundary_total_weights_u, i)
        ) / filtered_weights_determinant
      );
    }
  ); // same order as interior vertices
  const parametrized_interior_vertices_v = QUICK_VECTOR.vector_builder(
    interior_vertices_indices.length,
    function (i) {
      return (
        QUICK_MATRIX.determinant(
          QUICK_MATRIX.sub_column_for_vector(filtered_weights, boundary_total_weights_v, i)
        ) / filtered_weights_determinant
      );
    }
  ); // same order as interior vertices
  */

  const parametrized_interior_vertices = []; // same order as interior vertices
  for (let i = 0; i < parametrized_interior_vertices_u.length; i++) {
    parametrized_interior_vertices.push([
      parametrized_interior_vertices_u[i][0],
      parametrized_interior_vertices_v[i][0],
    ]);
  }

  const indices = interior_vertices_indices.concat(ordered_boundary_vertices_indices);
  const parametrized_vertices = parametrized_interior_vertices.concat(
    parametrized_ordered_boundary_vertices
  );

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

const auxetic_info = function (parametrized) {
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

  const max_rotation_angles = Math.PI / 5; // change later

  const max_radii = QUICK_VECTOR.apply_map(
    parametrized_vertices,
    function (parametrized_vertex, i) {
      const parametrized_vertex_neighbor_indices = neighbors(i, { edges: parametrized_edges });
      const lengths = [];
      for (const parametrized_vertex_neighbor_index of parametrized_vertex_neighbor_indices) {
        lengths.push(
          edge_length_2D([i, parametrized_vertex_neighbor_index], parametrized_vertices)
        );
      }
      const min_length = Math.min(...lengths);
      return min_length / 4;
    }
  );

  const radii = QUICK_VECTOR.apply_map(parametrized_edges, function (parametrized_edge, i) {
    const D =
      original_edge_lengths[i] - max_radii[parametrized_edge[0]] - max_radii[parametrized_edge[1]];
    const G = D - edge_length_difference[i] / 2;
    const a = 2 * (1 - Math.cos(max_rotation_angles));
    const b = 2 * G * (1 - Math.cos(max_rotation_angles));
    const c = G ** 2 - D ** 2;
    const r = (-b + Math.sqrt(b ** 2 - 4 * a * c)) / (2 * a);
    if (isNaN(r)) {
      return max_radii[parametrized_edge[0]];
    } else {
      return Math.max(
        Math.min(max_radii[parametrized_edge[0]], r),
        -max_radii[parametrized_edge[0]]
      );
    }
  });

  return {
    parametrized_vertices: parametrized_vertices,
    parametrized_edges: parametrized_edges,
    max_radii: max_radii,
    radii: radii,
    max_rotation_angles: max_rotation_angles,
  };
};

const auxetic_mesh = function (scene, models, auxetic_infoed) {
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

    /*
    const cloned_unit = unit.clone();
    cloned_unit.translateX(parametrized_vertex[0]);
    cloned_unit.translateY(parametrized_vertex[1]);
    cloned_unit.scale.set(max_radii[i], max_radii[i], boundary_lengthh / 40);
    scene.add(cloned_unit);
    */

    const neighbor_vertex_indices = neighbors(i, { edges: parametrized_edges });
    for (const neighbor_vertex_index of neighbor_vertex_indices) {
      const neighbor_vertex = parametrized_vertices[neighbor_vertex_index];
      const delta_x = neighbor_vertex[0] - parametrized_vertex[0];
      const delta_y = neighbor_vertex[1] - parametrized_vertex[1];
      const edge_length = edge_length_2D([i, neighbor_vertex_index], parametrized_vertices);
      const angle = Math.atan2(delta_y, delta_x) + max_rotation_angles;

      let wanted_edge_index;
      for (let j = 0; j < parametrized_edges.length; j++) {
        if (
          parametrized_edges[j].includes(i) === true &&
          parametrized_edges[j].includes(neighbor_vertex_index) === true
        ) {
          wanted_edge_index = j;
          break;
        }
      }

      const radius = radii[wanted_edge_index];

      if (radius >= 0.0001) {
        const cloned_rod_node = rod_mid.clone();
        cloned_rod_node.translateX(parametrized_vertex[0]);
        cloned_rod_node.translateY(parametrized_vertex[1]);
        cloned_rod_node.scale.set(radius, radius, boundary_lengthh / 10);
        cloned_rod_node.rotation.set(0, 0, angle - Math.PI / 2);
        scene.add(cloned_rod_node);

        const cloned_knob = knob.clone();
        cloned_knob.translateX(parametrized_vertex[0]);
        cloned_knob.translateY(parametrized_vertex[1]);
        const small_delta_x = radius * Math.cos(angle);
        const small_delta_y = radius * Math.sin(angle);
        cloned_knob.translateX(small_delta_x);
        cloned_knob.translateY(small_delta_y);
        cloned_knob.rotation.set(0, 0, angle - Math.PI / 2, "YXZ");
        cloned_knob.scale.set(radius, radius * 1.1, boundary_lengthh / 60);
        scene.add(cloned_knob);

        const average_x = (parametrized_vertex[0] + neighbor_vertex[0]) / 2;
        const average_y = (parametrized_vertex[1] + neighbor_vertex[1]) / 2;

        const base_x = parametrized_vertex[0] + small_delta_x + 0.55 * radius * Math.cos(angle);
        const base_y = parametrized_vertex[1] + small_delta_y + 0.55 * radius * Math.sin(angle);

        const rod_length = Math.sqrt((average_x - base_x) ** 2 + (average_y - base_y) ** 2) * 1.1;

        const cloned_rod_mid = rod_mid.clone();
        cloned_rod_mid.translateX(base_x);
        cloned_rod_mid.translateY(base_y);
        cloned_rod_mid.scale.set(radius / 2, rod_length, boundary_lengthh / 30);
        cloned_rod_mid.rotation.set(
          0,
          0,
          Math.atan2(average_y - base_y, average_x - base_x) - Math.PI / 2
        );
        scene.add(cloned_rod_mid);

        const cloned_rod_ext = rod_ext.clone();
        cloned_rod_ext.translateX(parametrized_vertex[0]);
        cloned_rod_ext.translateY(parametrized_vertex[1]);
        cloned_rod_ext.translateX(small_delta_x);
        cloned_rod_ext.translateY(small_delta_y);
        cloned_rod_ext.rotation.set(0, 0, angle - Math.PI / 2);
        cloned_rod_ext.scale.set(radius, radius * 1.1, boundary_lengthh / 30);
        scene.add(cloned_rod_ext);
      }
    }
  }
};

fetch("/public/untitled.obj").then(async function (response) {
  const test_obj = await response.text();
  const graphed = mesh_to_graph(test_obj);
  console.log("graphed");
  const boundaried = boundary(graphed);
  console.log("boundaried");
  const parametrized = parametrization(graphed, boundaried);
  console.log("parametrized");
  const auxetic_infoed = auxetic_info(parametrized);
  console.log("auxetic_infoed");

  const plotted = document.getElementById("plotted");
  const ctx = plotted.getContext("2d");
  const scale = 700;
  for (let j = 0; j < parametrized.parametrized_vertices.length; j++) {
    let radius;
    for (let i = 0; i < auxetic_infoed.parametrized_edges.length; i++) {
      if (auxetic_infoed.parametrized_edges[i].includes(j)) {
        radius = auxetic_infoed.radii[i];
      }
    }
    if (radius >= 0) {
      ctx.fillStyle = "black";
    } else {
      ctx.fillStyle = "grey";
    }
    ctx.beginPath();
    ctx.arc(
      parametrized.parametrized_vertices[j][0] * scale + 200,
      parametrized.parametrized_vertices[j][1] * scale + 100,
      Math.abs(radius) * scale,
      0,
      Math.PI * 2
    );
    ctx.fill();
  }
  for (const parametrized_edge of auxetic_infoed.parametrized_edges) {
    ctx.beginPath();
    ctx.moveTo(
      parametrized.parametrized_vertices[parametrized_edge[0]][0] * scale + 200,
      parametrized.parametrized_vertices[parametrized_edge[0]][1] * scale + 100
    );
    ctx.lineTo(
      parametrized.parametrized_vertices[parametrized_edge[1]][0] * scale + 200,
      parametrized.parametrized_vertices[parametrized_edge[1]][1] * scale + 100
    );
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = "red";
    ctx.stroke();
  }

  const models = {};
  {
    const loader = new GLTFLoader();
    loader.load("/public/unit.glb", function (gltf) {
      const mesh = gltf.scene;
      models.unit = mesh;
    });
    loader.load("/public/knob.glb", function (gltf) {
      const mesh = gltf.scene;
      models.knob = mesh;
    });
    loader.load("/public/rod-ext.glb", function (gltf) {
      const mesh = gltf.scene;
      models.rod_ext = mesh;
    });
    loader.load("/public/rod-mid.glb", function (gltf) {
      const mesh = gltf.scene;
      models.rod_mid = mesh;
    });
  }

  const original = document.getElementById("original");

  const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 999999999999999);
  camera.position.set(0.3, 0.38, 1);
  camera.rotation.set(0, 0, 0, "YXZ");

  const scene = new THREE.Scene();

  const renderer1 = new THREE.WebGLRenderer({
    canvas: original,
  });

  {
    const light = new THREE.PointLight(0xfff2e3, 30);
    light.position.set(2, 2, 5);
    scene.add(light);
  }

  {
    const sky = new Sky();
    sky.scale.setScalar(450000);

    const phi = Math.PI * 0.48;
    const theta = (Math.PI * 3) / 2;
    const sunPosition = new THREE.Vector3().setFromSphericalCoords(1, phi, theta);

    sky.material.uniforms.sunPosition.value = sunPosition;

    scene.add(sky);
  }

  {
    const loader = new OBJLoader();
    const response = loader.parse(test_obj);
    response.rotation.set(0, Math.PI / 4, 0, "YXZ");
    response.scale.set(-1, 1, 1);
    response.translateX(0.7);
    for (const child of response.children) {
      console.log(child);
      const material = new THREE.MeshPhongMaterial({ color: "grey" });
      child.material = material;
      child.material.side = THREE.DoubleSide;
    }
    scene.add(response);

    renderer1.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.zoom = Math.min(window.innerWidth / 1000, window.innerHeight / 1000);
    camera.updateProjectionMatrix();
    renderer1.render(scene, camera);

    scene.remove(response);
  }

  const canvas = document.getElementById("canvas");

  const renderer2 = new THREE.WebGLRenderer({
    canvas: canvas,
  });

  const waiter = function () {
    if (
      models.unit == undefined &&
      models.knob == undefined &&
      models.rod_ext == undefined &&
      models.rod_mid == undefined
    ) {
      setTimeout(waiter, 100);
    } else {
      auxetic_mesh(scene, models, auxetic_infoed);

      renderer2.setSize(window.innerWidth, window.innerHeight);
      camera.aspect = window.innerWidth / window.innerHeight;
      camera.zoom = Math.min(window.innerWidth / 1000, window.innerHeight / 1000);
      camera.updateProjectionMatrix();
      renderer2.render(scene, camera);

      const exporter = new GLTFExporter();
      exporter.parse(scene, function (gltf) {
        const blob = new Blob([JSON.stringify(gltf)], {
          type: "application/json",
        });
        const fileURL = URL.createObjectURL(blob);

        const link = document.createElement("a");
        link.href = fileURL;
        link.download = "auxetic-mesh.glb";
        link.innerHTML = "DOWNLOAD";
        document.body.appendChild(link);
      });
    }
  };
  waiter();
});
