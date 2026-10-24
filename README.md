<div align="center">

# 📐 Analytical Plate Solutions in MATLAB

**Navier and Lévy solutions for rectangular isotropic, orthotropic and laminated plates**

![MATLAB](https://img.shields.io/badge/MATLAB-0076A8?logo=mathworks&logoColor=white)
![Topic](https://img.shields.io/badge/topic-plate%20theory-blueviolet)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-archive-lightgrey)

[About](#-about) · [Contents](#-contents) · [Naming convention](#-naming-convention) · [Methods](#-methods) · [How to run](#-how-to-run) · [Contact](#-contact)

</div>

---

## 📌 About

This repository contains MATLAB codes with **analytical solutions for the bending of rectangular plates**, using the **Navier** and **Lévy** methods. The codes cover:

- **Single-layer plates**, isotropic and orthotropic
- **Laminated plates** with orthotropic layers
- **Sinusoidal** and **uniform** transverse loads

They can be used as **reference solutions**, for example to verify finite element models, and as **study material** on plate theory.

<!-- Add here the plate theory used (for example classical, first-order shear deformation) and the context of the codes (coursework, research) -->

## 📂 Contents

| File | Method | Plate | Load |
|---|---|---|---|
| [`Navier_solution.m`](Navier_solution.m) | Navier | General Navier solution routine | |
| [`Navier_single_iso_sin.m`](Navier_single_iso_sin.m) | Navier | Single layer, isotropic | Sinusoidal |
| [`Navier_single_iso_uni.m`](Navier_single_iso_uni.m) | Navier | Single layer, isotropic | Uniform |
| [`Navier_single_ort_sin.m`](Navier_single_ort_sin.m) | Navier | Single layer, orthotropic | Sinusoidal |
| [`Navier_single_ort_uni.m`](Navier_single_ort_uni.m) | Navier | Single layer, orthotropic | Uniform |
| [`Navier_laminado_ort_sin.m`](Navier_laminado_ort_sin.m) | Navier | Laminated, orthotropic layers | Sinusoidal |
| [`Navier_laminado_ort_uni.m`](Navier_laminado_ort_uni.m) | Navier | Laminated, orthotropic layers | Uniform |
| [`Levy_single_ort_uni.m`](Levy_single_ort_uni.m) | Lévy | Single layer, orthotropic | Uniform |

## 🏷️ Naming convention

File names follow the pattern `Method_plate_material_load.m`:

| Part | Meaning |
|---|---|
| `Navier` / `Levy` | Analytical method |
| `single` / `laminado` | Single-layer plate / laminated plate (*laminado* is Portuguese for laminated) |
| `iso` / `ort` | Isotropic / orthotropic material |
| `sin` / `uni` | Sinusoidal / uniform load |

For example, `Navier_laminado_ort_uni.m` is the Navier solution for a laminated orthotropic plate under uniform load.

## 🧠 Methods

- **Navier solution**: for rectangular plates **simply supported on all four edges**. The load and the displacements are expanded in double Fourier series, which turns the governing equations into algebraic equations for each term of the series.
- **Lévy solution**: for rectangular plates with **two opposite edges simply supported**, while the other two edges can have different boundary conditions (free, clamped or simply supported). The solution uses a single Fourier series, and the remaining ordinary differential equation is solved analytically.

<!--
## 🖼️ Results
Add plots of the deflection and stresses here. Save the images in an images/ folder, for example:

| Isotropic plate | Orthotropic plate | Laminated plate |
|:---:|:---:|:---:|
| ![Iso](images/iso.png) | ![Ort](images/ort.png) | ![Lam](images/laminated.png) |
-->

## 🚀 How to run

**Requirements:** MATLAB.

```matlab
% 1. Clone or download this repository and open the folder in MATLAB
% 2. Run any script by typing its name in the Command Window:
Navier_single_iso_uni
```

You can also open a script in the MATLAB editor and click **Run**.

<!-- Describe the inputs (plate dimensions, material properties, load, number of series terms) and the outputs here -->

## 📚 Further reading

- Reddy, J. N. (2004). *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed. CRC Press.
- Timoshenko, S. and Woinowsky-Krieger, S. (1959). *Theory of Plates and Shells*, 2nd ed. McGraw-Hill.

## 🤝 How to contribute

Suggestions and fixes are welcome. Open an [issue](../../issues) or send a **Pull Request**.

## 🔗 Related projects

- [fem-python-lab](https://github.com/RenanMP12/fem-python-lab): FEM and numerical methods scripts in Python
- [AE-245](https://github.com/RenanMP12/AE-245): Finite Element Method coursework at ITA, in Python

## 📜 License

Distributed under the **MIT** license. See the [`LICENSE`](LICENSE) file for details.

## 📬 Contact

**Renan Portela**, PhD in Mechanical Engineering (University of Waterloo, Canada) · Master's in Mechanical Engineering (UNICAMP)

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0A66C2?logo=linkedin&logoColor=white)](https://www.linkedin.com/in/YOUR-PROFILE)
[![Email](https://img.shields.io/badge/E--mail-D14836?logo=gmail&logoColor=white)](mailto:YOUR-EMAIL@example.com)

---

<div align="center">

⭐ If this repository helped you, leave a star!

</div>
