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
    coefficients = 'coef_half coef_half'
    powers = '3 3'
    base_name = lps
    outputs = all
    initial_porosity = 0.1
  [../]
  [./coef_half]
    type = ParsedMaterial
    args = temp
    f_name = coef_half
    function = '4e-20 * exp(-28500 / 1.987 / 1200)'
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

  end_time = 1e-1
  dt = 1e-2
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
  [./lps_eff_creep_strain]
    type = ElementAverageValue
    variable = lps_effective_creep_strain
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = porosity
    execute_on = 'TIMESTEP_END initial'
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
  perf_graph = true
[]
