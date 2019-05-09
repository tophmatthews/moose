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
  volumetric_locking_correction = false
[]

[AuxVariables]
  [./temp]
    initial_condition = 800
  [../]
[]

[Modules/TensorMechanics/Master/All]
  strain = FINITE
  add_variables = true
  generate_output = 'strain_xx strain_yy strain_xy hydrostatic_stress vonmises_stress'
  use_automatic_differentiation = true
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
    type = ADComputeMultipleInelasticStress
    inelastic_models = lps
  [../]

  [./lps]
    type = ADLPSViscoPlasticityStressUpdate
    activation_energy = 28500
    n_exponent = 3
    coefficient = 8e-20
    gas_constant = 1.987
    temperature = temp
    base_name = lps
    outputs = all
    initial_porosity = 1e-1
  [../]
[]

[BCs]
  [./no_disp_x]
    type = ADPresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./no_disp_y]
    type = ADPresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./pull_disp_y]
    type = ADFunctionPresetBC
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
  line_search = 'none'

  end_time = 2
  timestep_tolerance = 1e-6

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
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
  [../]
  [./num_nonlin]
    type = NumNonlinearIterations
  [../]
  [./lps_eff_creep_strain]
    type = ElementAverageValue
    variable = lps_effective_creep_strain
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = porosity
    execute_on = 'TIMESTEP_END initial'
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
  execute_on = 'TIMESTEP_END'
  exodus = true
  print_linear_residuals = false
  perf_graph = true
  csv = true
[]
