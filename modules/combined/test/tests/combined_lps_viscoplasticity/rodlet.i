[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 0.002
  ymax = 0.02
  second_order = true
[]

[Problem]
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  order = SECOND
[]

[Variables]
  [./temp]
    initial_condition = 600
  [../]
[]

[Kernels]
  [./heat_conduction]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = 20.0
  [../]
  [./heat_dt]
    type = ADHeatConductionTimeDerivative
    variable = temp
    specific_heat = 1.0
  [../]
  [./heat_source]
    type = ADBodyForce
    variable = temp
    function = temp_ramp
  [../]
[]

[Functions]
  [./temp_ramp]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 1e10 1e10'
  [../]
[]

[Modules/TensorMechanics/Master/All]
  strain = FINITE
  add_variables = true
  generate_output = 'strain_xx strain_yy strain_xy hydrostatic_stress vonmises_stress'
  eigenstrain_names = 'thermal_expansion'
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
  [./thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    eigenstrain_name = thermal_expansion
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    temperature = temp
  [../]
  [./density]
    type = Density
    density = 1.0e5
  [../]

  [./lps]
    type = LPSViscoPlasticityStressUpdate
    coefficients = coef
    powers = 3
    base_name = lps
    outputs = all
    initial_porosity = 0.1
    relative_tolerance = 1e-11
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

  [./temp]
    type = PresetBC
    variable = temp
    value = 600
    boundary = 'bottom right top'
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
  line_search = 'none'
  # l_tol = 1e-03
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm       lu'

  end_time = 2
  timestep_tolerance = 1e-6
  nl_max_its = 20

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    timestep_limiting_postprocessor = creep_timestep
  [../]
[]

[Postprocessors]
  [./max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
  [../]
  [./max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
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
  [./porosity_max]
    type = ElementExtremeValue
    variable = porosity
  [../]
  [./porosity_min]
    type = ElementExtremeValue
    variable = porosity
    value_type = min
  [../]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  [../]
  [./creep_timestep]
    type = MaterialTimeStepPostprocessor
  [../]
  [./temp_avg]
    type = ElementAverageValue
    variable = temp
  [../]
  [./temp_max]
    type = ElementExtremeValue
    variable = temp
  [../]
  [./temp_min]
    type = ElementExtremeValue
    variable = temp
    value_type = min
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
