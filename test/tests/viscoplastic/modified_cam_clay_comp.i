[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = 0
  zmax = 0.1
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mech_x]
    type = VPStressDivergence
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = VPStressDivergence
    variable = disp_y
    component = 1
  [../]
  [./mech_z]
    type = VPStressDivergence
    variable = disp_z
    component = 2
  [../]
[]

[AuxVariables]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dev_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_vol_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_eqv_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./yield]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pressure_aux]
    type = VPPressureAux
    variable = pressure
    execute_on = 'TIMESTEP_END'
  [../]
  [./dev_stress_aux]
    type = VPVonMisesStressAux
    variable = dev_stress
    execute_on = 'TIMESTEP_END'
  [../]
  [./plastic_vol_strain_aux]
    type = VPVolStrainAux
    variable = plastic_vol_strain
    strain_type = 'plastic'
    execute_on = 'TIMESTEP_END'
  [../]
  [./plastic_eqv_strain_aux]
    type = VPEqvStrainAux
    variable = plastic_eqv_strain
    strain_type = 'plastic'
    execute_on = 'TIMESTEP_END'
  [../]
  [./yield_aux]
    type = MaterialRealAux
    variable = yield
    property = yield_function
    execute_on = 'TIMESTEP_END'
  [../]
[]

[BCs]
  [./x_compression]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 'left right'
    function = '-1.0e-05*x*t'
  [../]
  [./y_extension]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 'bottom top'
    function = '1.0e-05*y*t'
  [../]
  [./no_z]
    type = PresetBC
    variable = disp_z
    boundary = 'front back'
    value = 0.0
  [../]
[]

[Materials]
  [./elastic_mat]
    type = VPMechMaterial
    displacements = 'disp_x disp_y disp_z'
    bulk_modulus = 2.0e+03
    shear_modulus = 2.0e+03
    initial_stress = '-7.0 -7.0 -7.0'
    viscoplastic_model = 'mod_cam_clay_mat'
  [../]
  [./mod_cam_clay_mat]
    type = VPMCamClay
    friction_angle = 30.0
    critical_pressure = 1.0e+01
    plastic_viscosity = 2.0e+04
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -snes_atol'
    petsc_options_value = 'hypre boomeramg 0.7 4 5 25 HMIS ext+i 2 0.3 1.0e-08'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  start_time = 0.0
  end_time = 100.0
  dt = 5.0
[]

[Outputs]
  execute_on = 'TIMESTEP_END'
  print_linear_residuals = true
  perf_graph = true
  exodus = true
  # csv = true
[]
