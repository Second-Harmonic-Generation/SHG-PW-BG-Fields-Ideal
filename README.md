## 🧰 How to Use This Template    

Click the green **"Use this template"** button at the top of the page, then choose **"Create a new repository"**.   

This will create your own copy of this project, which you can modify freely — no need to fork!   

 
<p align="center">
  <img src="./images/SHG-banner.png" alt="SHG Logo">
</p>


<h1 align="center">SHG-PW-BG-Field-Ideal</h1>

<div align="center">

| **Term** | **Definition** |
|----------|----------------|
| **SHG** | Second Harmonic Generation |
| **PW** | Pulsed Wave |
| **BG** | Bessel-Gaussian |
</div>

&nbsp;

<div align="center">

Article title:       
**Pulsed Bessel–Gauss beams: a depleted wave model for type II second-harmonic generation**
</div>

&nbsp;

---

***Table of Contents***

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#1-about-this-repository"><i><b>1. About this repository</b></i></a>
</div>
&nbsp;

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#2-getting-started"><i><b>2. Getting Started</b></i></a>
</div>
&nbsp;

<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#3-how-to-cite-us"><i><b>3. How to Cite Us</b></i></a>
</div>
&nbsp;


<div>
  &nbsp;&nbsp;&nbsp;&nbsp;<a href="#4-contact-information"><i><b>4. Contact Information</b></i></a>
</div>
&nbsp;

---    

## 1. About this repository

This repository contains the **Toolkit for Modeling of 3D Temperature Distribution in KTP Crystal: Continuous-Wave Gaussian Second Harmonic Generation**, an open-source toolkit for modeling the thermal dynamics that underpin continuous-wave second-harmonic generation (CW SHG), using KTP as a case study.

### Toolkit Overview

The toolkit provides comprehensive modules for geometry and material definitions of KTP crystals, boundary and cooling models with various heat transfer mechanisms, and transient and steady-state finite-difference solvers for temperature field computation.

The toolkit supports parameterized scenario sweeps including temperature-dependent versus constant thermal conductivity, convection with and without radiation boundary conditions, and heat-transfer coefficients spanning 6.5–2.0×10⁴ W·m⁻²·K⁻¹. It features compiled Fortran kernels with built-in benchmark reporting, reproducible pipelines with versioned code repository, and exportable datasets with spatiotemporal temperature fields. The toolkit generates both radial and axial temperature profiles for comprehensive analysis.

The implementation has been validated by reproducing temperature distributions and trends for KTP under Gaussian CW pumping, including the effects of temperature-dependent conductivity and boundary conditions. This toolkit was used to solve the thermal modeling problem described in the research article **"Temperature Distribution in a Gaussian End-Pumped Nonlinear KTP Crystal: the Temperature Dependence of Thermal Conductivity and Radiation Boundary Condition"**.  


```
Folder PATH listing
+---citation                      <-- Contains citation materials and research papers
│       1_Heat-Equation_Continu…  <-- Heat equation analytical paper
│       2_Heat-Equation_Continu…  <-- Heat equation continuous wave paper
│       3_Heat-Equation_Pulsed-…  <-- Heat equation pulsed wave paper
│       4_Phase-Mismatch_Pulsed…  <-- Phase mismatch pulsed wave paper
│       5_Ideal_Continuous-Wave…  <-- Ideal continuous wave paper
│       6_Ideal_Pulsed-Wave_Be…   <-- Ideal pulsed wave Bessel paper
│       7_Coupled_Continuous-Wa…  <-- Coupled continuous wave paper
│       README.md                 <-- Citation guidelines and information
│
+---images                        <-- Contains project images and logos
│       SHG-banner.png            <-- SHG project banner
│
+---results                       <-- Numerical simulation results and benchmark data
│       E08_f_4000_Np_1_tp_50_…   <-- Best iteration results data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component r-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component t-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component z-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component r-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component t-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component z-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component r-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component t-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Electric field component z-direction data
│       E_08_f_4000_Np_1_tp_50…   <-- Scale maximum data
│       E_08_f_4000_Np_1_tp_50…   <-- Phase pick data for Psi2
│       E_08_f_4000_Np_1_tp_50…   <-- Phase pick data for Psi3
│
+---src                           <-- Toolkit source code and documentation
│       Code_SHG-PW-BG-Field-…    <-- Fortran finite difference solver (main toolkit)
│
│       Article_SHG-PW-BG-Fiel…   <-- Research paper PDF (problem solved by toolkit)
│       CITATION.cff              <-- Citation metadata file
│       LICENSE                   <-- Project license information
│       README.md                 <-- Toolkit overview and documentation
│
```

## 2. Getting Started

### 2.1. Prerequisites
- **Fortran Compiler** (gfortran, Intel Fortran, or similar)
- **Text Editor** (VS Code, Cursor, or any Fortran-capable editor)
- **PDF Reader** (for accessing research papers)
- **Git** (for cloning the repository)
- **Plotting Software** (Gnuplot, Python matplotlib, or similar for visualizing results)

### 2.2. Quick Start

1. **Clone the Repository**
   ```bash
   git clone https://github.com/Second-Harmonic-Generation/SHG-PW-BG-Fields-Ideal.git
   cd SHG-PW-BG-Fields-Ideal
   ```

2. **Explore the Research Papers**
   - Navigate to the `citation/` folder
   - Read the main research paper: `Article_SHG-PW-BG-Field-Ideal.pdf`
   - Review additional papers for comprehensive understanding

3. **Compile and Run the Simulation**
   ```bash
   cd src
   gfortran -o shg_simulation Code_SHG-PW-BG-Field-Ideal.f90
   ./shg_simulation
   ```

4. **Analyze Results**
   - Check the `results/` folder for output files
   - Examine `.plt` files for electric field component data
   - Use plotting software to visualize the simulation results:
     ```bash
     # Example with Gnuplot
     gnuplot -e "plot 'results/E_08_f_4000_Np_1_tp_50_Elec12_r.plt' with lines"
     ```

5. **Explore Different Parameters** (Optional)
   - Edit the Fortran source code to modify simulation parameters
   - Recompile and run to explore different scenarios
   - Compare results with published findings in the research papers



## 3. How to Cite Us
Please refer to the [**citation**](./citation/) folder for accurate citations. It contains essential guidelines for accurate referencing, ensuring accurate acknowledgement of our work.


  
## 4. Contact Information

For questions not addressed in the resources above, please connect with [Mostafa Rezaee](https://www.linkedin.com/in/mostafa-rezaee/) on LinkedIn for personalized assistance.
