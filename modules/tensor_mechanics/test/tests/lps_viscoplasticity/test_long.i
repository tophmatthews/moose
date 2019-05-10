[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
  nz = 4
  xmax = 0.002
  ymax = 0.002
  zmax = 0.002
  second_order = true
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules/TensorMechanics/Master/All]
  strain = FINITE
  add_variables = true
  generate_output = 'strain_xx strain_yy strain_xy hydrostatic_stress vonmises_stress'
[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1e-4'
  [../]
[]


[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e10
    poissons_ratio = 0.3
  [../]
  [./strain]
    type = ComputeMultipleInelasticStress
    inelastic_models = lps
  [../]

  [./lps]
    type = LPSViscoPlasticityStressUpdate
    coefficients = coef
    powers = 3
    base_name = lps
    outputs = all
    initial_porosity = 0.1
    relative_tolerance = 1e-13
  [../]
  [./coef]
    type = ParsedMaterial
    f_name = coef
    function = '8e-20 * exp(-28500 / 1.987 / 1200)'
  [../]
[]

[BCs]
  [./no_disp_x]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./no_disp_y]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./pull_disp_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = pull
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  end_time = 2
  timestep_tolerance = 1e-6

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    timestep_limiting_postprocessor = creep_timestep
  [../]
[]

[Postprocessors]
  [./disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  [../]
  [./disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./max_hydro]
    type = ElementExtremeValue
    variable = hydrostatic_stress
  [../]
  [./min_hydro]
    type = ElementExtremeValue
    variable = hydrostatic_stress
    value_type = min
  [../]
  [./avg_hydro]
    type = ElementAverageValue
    variable = hydrostatic_stress
  [../]
  [./max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
  [../]
  [./min_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = min
  [../]
  [./avg_vonmises]
    type = ElementAverageValue
    variable = vonmises_stress
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./num_lin]
    type = NumLinearIterations
    outputs = console
  [../]
  [./num_nonlin]
    type = NumNonlinearIterations
    outputs = console
  [../]
  [./total_lin]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin
    outputs = console
  [../]
  [./total_nonlin]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin
    outputs = console
  [../]
  [./lps_eff_creep_strain]
    type = ElementAverageValue
    variable = lps_effective_creep_strain
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = porosity
  [../]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  [../]
  [./creep_timestep]
    type = MaterialTimeStepPostprocessor
  [../]
[]

[Outputs]
  perf_graph = true
  [./out]
    type = CSV
    sync_only = true
    sync_times = '0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0'
  [../]
[]
