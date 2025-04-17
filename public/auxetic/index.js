import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js";
import { OBJLoader } from "three/addons/loaders/OBJLoader.js";
import { GLTFExporter } from "three/addons/exporters/GLTFExporter.js";

import {
  mesh_to_graph,
  boundary,
  parametrization,
  auxetic_info,
  auxetic_mesh,
  folded_mesh,
  morph_intestine,
  morph_forces,
} from "/public/auxetic/auxetics.js";

document.getElementById("3d-obj").addEventListener("change", function () {
  const file = document.getElementById("3d-obj").files[0];

  const reader = new FileReader();
  reader.addEventListener("load", function () {
    const width = parseInt(document.getElementById("width").value);
    const height = parseInt(document.getElementById("height").value);
    const starting_vertex = parseInt(document.getElementById("starting").value);
    print(
      reader.result,
      width,
      height,
      starting_vertex,
      600,
      {
        x: 0,
        y: 0,
        z: 1,
      },
      400
    );
  });

  reader.readAsText(file);
});

/*
fetch("/public/auxetic/heights.obj").then(async function (response) {
  const test_obj = await response.text();
  print(
    test_obj,
    5,
    5,
    1,
    200,
    {
      x: 0,
      y: 0.5,
      z: 4,
    },
    500
  );
});
*/

/*
fetch("/public/auxetic/dome.obj").then(async function (response) {
  const test_obj = await response.text();
  print(test_obj, 3, 3, 1, 200, {
    x: 0,
    y: 0,
    z: 4,
  }, 400);
});
*/

fetch("/public/auxetic/untitled.obj").then(async function (response) {
  const test_obj = await response.text();
  print(
    test_obj,
    10,
    6,
    1,
    600,
    {
      x: 0,
      y: 0,
      z: 1,
    },
    400
  );
});

const print_plot = function (boundaried, parametrized, auxetic_infoed, scale, second) {
  document.getElementById("boundary-size").innerHTML = "Boundary size: " + boundaried.ordered_boundary_vertices_indices.length;

  const plotted = document.getElementById("plotted");
  const ctx = plotted.getContext("2d");

  const highest_radius = Math.max(
    Math.abs(Math.max(...auxetic_infoed.parametrized_edge_length_differences)),
    Math.abs(Math.min(...auxetic_infoed.parametrized_edge_length_differences))
  );
  for (let i = 0; i < auxetic_infoed.parametrized_edges.length; i++) {
    const parametrized_edge = auxetic_infoed.parametrized_edges[i];
    let radius = auxetic_infoed.parametrized_edge_length_differences[i]; // prevent name repeat
    if (radius >= 0) {
      ctx.strokeStyle =
        "rgba(0, 0, " + String((Math.abs(radius) / highest_radius) * 255) + ", " + String(Math.abs(radius) / highest_radius + 0.1) + ")";
      ctx.lineWidth = (Math.abs(radius) / highest_radius) * 5 + 0.1;
    } else {
      ctx.strokeStyle =
        "rgba(" + String((Math.abs(radius) / highest_radius) * 255) + ", 0, 0, " + String(Math.abs(radius) / highest_radius + 0.1) + ")";
      ctx.lineWidth = (Math.abs(radius) / highest_radius) * 5 + 0.1;
    }
    ctx.beginPath();
    ctx.moveTo(
      parametrized.parametrized_vertices[parametrized_edge[0]][0] * scale,
      plotted.height - parametrized.parametrized_vertices[parametrized_edge[0]][1] * scale
    );
    ctx.lineTo(
      parametrized.parametrized_vertices[parametrized_edge[1]][0] * scale,
      plotted.height - parametrized.parametrized_vertices[parametrized_edge[1]][1] * scale
    );
    ctx.closePath();
    ctx.stroke();
  }

  for (const parametrized_edge of auxetic_infoed.parametrized_edges) {
    ctx.beginPath();
    ctx.moveTo(
      parametrized.parametrized_vertices[parametrized_edge[0]][0] * scale,
      plotted.height - (parametrized.parametrized_vertices[parametrized_edge[0]][1] * scale + second)
    );
    ctx.lineTo(
      parametrized.parametrized_vertices[parametrized_edge[1]][0] * scale,
      plotted.height - (parametrized.parametrized_vertices[parametrized_edge[1]][1] * scale + second)
    );
    ctx.lineWidth = 0.5;
    ctx.strokeStyle = "red";
    ctx.closePath();
    ctx.stroke();
  }
  for (let j = 0; j < parametrized.parametrized_vertices.length; j++) {
    let angle = auxetic_infoed.max_rotation_angles[j];
    let radius = auxetic_infoed.max_radii[j];

    if (angle >= 0) {
      ctx.fillStyle = "blue";
    } else {
      ctx.fillStyle = "grey";
    }
    ctx.beginPath();
    ctx.moveTo(parametrized.parametrized_vertices[j][0] * scale, plotted.height - (parametrized.parametrized_vertices[j][1] * scale + second));
    ctx.arc(
      parametrized.parametrized_vertices[j][0] * scale,
      plotted.height - (parametrized.parametrized_vertices[j][1] * scale + second),
      Math.abs(radius) * scale,
      0,
      angle
    );
    ctx.closePath();
    ctx.fill();
    ctx.beginPath();
    ctx.moveTo(parametrized.parametrized_vertices[j][0] * scale, plotted.height - (parametrized.parametrized_vertices[j][1] * scale + second));
    ctx.arc(
      parametrized.parametrized_vertices[j][0] * scale,
      plotted.height - (parametrized.parametrized_vertices[j][1] * scale + second),
      Math.abs(radius) * scale,
      0,
      2 * Math.PI
    );
    ctx.closePath();
    ctx.strokeStyle = "black";
    ctx.stroke();
  }
};

const print_original = function (text, cam_pos) {
  const original = document.getElementById("original");

  const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 999999999999999);
  camera.position.set(cam_pos.x, cam_pos.y, cam_pos.z);
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
    const loader = new OBJLoader();
    const response = loader.parse(text);
    for (const child of response.children) {
      const material = new THREE.MeshPhongMaterial({ color: "grey" });
      child.material = material;
      child.material.side = THREE.DoubleSide;
    }
    scene.add(response);

    const body_container = document.getElementById("body-container");

    renderer1.setSize(body_container.offsetWidth, body_container.offsetWidth);
    camera.aspect = body_container.offsetWidth / body_container.offsetWidth;
    camera.zoom = Math.min(body_container.offsetWidth / 1000, body_container.offsetWidth / 1000);
    camera.updateProjectionMatrix();

    let rot = 0;

    const perpetual = function () {
      rot += 0.01;
      response.rotation.set(0, rot, 0, "YXZ");

      renderer1.render(scene, camera);

      requestAnimationFrame(perpetual);
    };
    perpetual();
  }
};

const print_generated = function (models, auxetic_infoed) {
  const canvas = document.getElementById("canvas");

  const camera_position = {
    x: 0.5,
    y: 0.3,
    z: 0.5,
  };

  const camera = new THREE.PerspectiveCamera(40, 1, 0.0001, 999999999999999);
  camera.position.set(camera_position.x, camera_position.y, camera_position.z);
  camera.rotation.set(0, 0, 0, "YXZ");

  const scene = new THREE.Scene();

  const renderer = new THREE.WebGLRenderer({
    canvas: canvas,
  });

  {
    const light = new THREE.PointLight(0xfff2e3, 30);
    light.position.set(2, 2, 5);
    scene.add(light);
  }

  const body_container = document.getElementById("body-container");

  renderer.setSize(body_container.offsetWidth, body_container.offsetWidth);
  camera.aspect = body_container.offsetWidth / body_container.offsetWidth;
  camera.zoom = Math.min(body_container.offsetWidth / 1000, body_container.offsetWidth / 1000);
  camera.updateProjectionMatrix();
  renderer.render(scene, camera);

  const waiter = function () {
    if (models.unit === undefined && models.knob === undefined && models.rod_ext === undefined && models.rod_mid === undefined) {
      setTimeout(waiter, 300);
    } else {
      auxetic_mesh(scene, models, auxetic_infoed, 0.5);

      renderer.setSize(body_container.offsetWidth, body_container.offsetWidth);
      camera.aspect = body_container.offsetWidth / body_container.offsetWidth;
      camera.zoom = Math.min(body_container.offsetWidth / 1000, body_container.offsetWidth / 1000);
      camera.updateProjectionMatrix();
      renderer.render(scene, camera);

      const exporter = new GLTFExporter();
      exporter.parse(scene, function (gltf) {
        const blob = new Blob([JSON.stringify(gltf)], {
          type: "application/json",
        });
        const fileURL = URL.createObjectURL(blob);

        const link = document.createElement("a");
        link.href = fileURL;
        link.download = "auxetic-mesh.glb";
        link.innerHTML = "Download generated mesh";
        document.getElementById("download-div").appendChild(link);
      });

      const active_keys = {
        w: false,
        a: false,
        s: false,
        d: false,
        q: false,
        e: false,
      };

      window.addEventListener("keydown", function (e) {
        switch (e.key) {
          case "w":
            active_keys.w = true;
            break;
          case "a":
            active_keys.a = true;
            break;
          case "s":
            active_keys.s = true;
            break;
          case "d":
            active_keys.d = true;
            break;
          case "q":
            active_keys.q = true;
            break;
          case "e":
            active_keys.e = true;
            break;
        }
      });
      window.addEventListener("keyup", function (e) {
        switch (e.key) {
          case "w":
            active_keys.w = false;
            break;
          case "a":
            active_keys.a = false;
            break;
          case "s":
            active_keys.s = false;
            break;
          case "d":
            active_keys.d = false;
            break;
          case "q":
            active_keys.q = false;
            break;
          case "e":
            active_keys.e = false;
            break;
        }
      });

      const perpetual = function () {
        if (active_keys.q === true) {
          camera_position.z *= 0.95;
        }
        if (active_keys.e === true) {
          camera_position.z /= 0.95;
        }
        if (active_keys.a === true) {
          camera_position.x -= 0.01 * camera_position.z;
        }
        if (active_keys.d === true) {
          camera_position.x += 0.01 * camera_position.z;
        }
        if (active_keys.w === true) {
          camera_position.y += 0.01 * camera_position.z;
        }
        if (active_keys.s === true) {
          camera_position.y -= 0.01 * camera_position.z;
        }
        if (
          active_keys.w === true ||
          active_keys.a === true ||
          active_keys.s === true ||
          active_keys.d === true ||
          active_keys.q === true ||
          active_keys.e === true
        ) {
          camera.position.set(camera_position.x, camera_position.y, camera_position.z);
          camera.updateProjectionMatrix();
          renderer.render(scene, camera);
        }
        requestAnimationFrame(perpetual);
      };
      perpetual();
    }
  };
  waiter();
};

const print = function (text, width, height, starting, scale, cam_pos, second) {
  const graphed = mesh_to_graph(text);
  const boundaried = boundary(graphed, width, height, starting);
  const parametrized = parametrization(graphed, boundaried, function (edge_length) {
    return 1 / edge_length ** 2;
  });
  const auxetic_infoed = auxetic_info(parametrized, 0.4);

  print_original(text, cam_pos);
  print_plot(boundaried, parametrized, auxetic_infoed, scale, second);

  const models = {};
  {
    const loader = new GLTFLoader();
    loader.load("/public/auxetic/unit.glb", function (gltf) {
      const mesh = gltf.scene;
      models.unit = mesh;
    });
    loader.load("/public/auxetic/knob.glb", function (gltf) {
      const mesh = gltf.scene;
      models.knob = mesh;
    });
    loader.load("/public/auxetic/rod-ext.glb", function (gltf) {
      const mesh = gltf.scene;
      models.rod_ext = mesh;
    });
    loader.load("/public/auxetic/rod-mid.glb", function (gltf) {
      const mesh = gltf.scene;
      models.rod_mid = mesh;
    });
  }

  print_generated(models, auxetic_infoed);
  print_folded(models, parametrized, auxetic_infoed, graphed);
};

const print_folded = function (models, parametrized, auxetic_infoed, graphed) {
  const canvas = document.getElementById("folded");

  const camera_position = {
    x: 0.5,
    y: 0.3,
    z: 0.5,
  };

  const camera = new THREE.PerspectiveCamera(40, 1, 0.0001, 999999999999999);
  camera.position.set(camera_position.x, camera_position.y, camera_position.z);
  camera.rotation.set(0, 0, 0, "YXZ");

  const scene = new THREE.Scene();

  const renderer = new THREE.WebGLRenderer({
    canvas: canvas,
  });

  {
    const light = new THREE.PointLight(0xfff2e3, 30);
    light.position.set(2, 2, 5);
    scene.add(light);
  }

  const body_container = document.getElementById("body-container");

  renderer.setSize(body_container.offsetWidth, body_container.offsetWidth);
  camera.aspect = body_container.offsetWidth / body_container.offsetWidth;
  camera.zoom = Math.min(body_container.offsetWidth / 1000, body_container.offsetWidth / 1000);
  camera.updateProjectionMatrix();
  renderer.render(scene, camera);

  const waiter = function () {
    if (models.unit === undefined && models.knob === undefined && models.rod_ext === undefined && models.rod_mid === undefined) {
      setTimeout(waiter, 300);
    } else {
      const folded_models = folded_mesh(scene, models, parametrized, auxetic_infoed, 0.5);

      // gravity and rotations
      const external_forces = [];
      const unit_rotations = [];
      for (let i = 0; i < graphed.vertices.length; i++) {
        external_forces.push([0, -1, 0]);
        unit_rotations.push(Math.PI / 6);
      }

      try {
        const morph_forced = morph_forces(graphed, parametrized, auxetic_infoed, external_forces, unit_rotations);
      } catch {}

      renderer.setSize(body_container.offsetWidth, body_container.offsetWidth);
      camera.aspect = body_container.offsetWidth / body_container.offsetWidth;
      camera.zoom = Math.min(body_container.offsetWidth / 1000, body_container.offsetWidth / 1000);
      camera.updateProjectionMatrix();
      renderer.render(scene, camera);

      const exporter = new GLTFExporter();
      exporter.parse(scene, function (gltf) {
        const blob = new Blob([JSON.stringify(gltf)], {
          type: "application/json",
        });
        const fileURL = URL.createObjectURL(blob);

        const link = document.createElement("a");
        link.href = fileURL;
        link.download = "auxetic-mesh-folded.glb";
        link.innerHTML = " Download folded mesh";
        document.getElementById("download-div").appendChild(link);
      });

      const active_keys = {
        w: false,
        a: false,
        s: false,
        d: false,
        q: false,
        e: false,
      };

      window.addEventListener("keydown", function (e) {
        switch (e.key) {
          case "w":
            active_keys.w = true;
            break;
          case "a":
            active_keys.a = true;
            break;
          case "s":
            active_keys.s = true;
            break;
          case "d":
            active_keys.d = true;
            break;
          case "q":
            active_keys.q = true;
            break;
          case "e":
            active_keys.e = true;
            break;
        }
      });
      window.addEventListener("keyup", function (e) {
        switch (e.key) {
          case "w":
            active_keys.w = false;
            break;
          case "a":
            active_keys.a = false;
            break;
          case "s":
            active_keys.s = false;
            break;
          case "d":
            active_keys.d = false;
            break;
          case "q":
            active_keys.q = false;
            break;
          case "e":
            active_keys.e = false;
            break;
        }
      });

      let old;
      let current = performance.now();
      const perpetual = function () {
        old = current;
        current = performance.now();
        let delta = Math.min(current - old, 20);
        if (active_keys.q === true) {
          camera_position.z *= 0.95;
        }
        if (active_keys.e === true) {
          camera_position.z /= 0.95;
        }
        if (active_keys.a === true) {
          camera_position.x -= 0.01 * camera_position.z;
        }
        if (active_keys.d === true) {
          camera_position.x += 0.01 * camera_position.z;
        }
        if (active_keys.w === true) {
          camera_position.y += 0.01 * camera_position.z;
        }
        if (active_keys.s === true) {
          camera_position.y -= 0.01 * camera_position.z;
        }

        morph_intestine(parametrized, auxetic_infoed, folded_models, current / 1000, delta / 1000, 0.5);

        camera.position.set(camera_position.x, camera_position.y, camera_position.z);
        camera.updateProjectionMatrix();
        renderer.render(scene, camera);

        requestAnimationFrame(perpetual);
      };
      perpetual();
    }
  };
  waiter();
};
