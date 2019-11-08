E_block = 1e7
E_plank = 1e7
elem = QUAD4
order = FIRST
name = 'small'

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./plank]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -0.3
    xmax = 0.3
    ymin = -10
    ymax = 10
    nx = 5
    ny = 67
    elem_type = ${elem}
  [../]
  [./plank_sidesets]
    type = RenameBoundaryGenerator
    input = plank
    old_boundary_id = '0 1 2 3'
    new_boundary_name = 'plank_bottom plank_right plank_top plank_left'
  [../]
  [./plank_id]
    type = SubdomainIDGenerator
    input = plank_sidesets
    subdomain_id = 1
  [../]

  [./block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.31
    xmax = 0.91
    ymin = 7.7
    ymax = 8.5
    nx = 5
    ny = 4
    elem_type = ${elem}
  [../]
  [./block_id]
    type = SubdomainIDGenerator
    input = block
    subdomain_id = 2
  [../]

  [./combined]
    type = MeshCollectionGenerator
    inputs = 'plank_id block_id'
  [../]
  [./block_rename]
    type = RenameBlockGenerator
    input = combined
    old_block_id = '1 2'
    new_block_name = 'plank block'
  [../]
  [./block_sidesets]
    type = SideSetsFromPointsGenerator
    input = block_rename
    points = '0.6  7.7  0
              0.91 8.0  0
              0.6  8.5 0
              0.31 8.0  0'
    new_boundary = 'block_bottom block_right block_top block_left'
  [../]

  [./slave]
    input = block_sidesets
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'block_left'
    new_block_id = '30'
    new_block_name = 'frictionless_slave_subdomain'
  [../]
  [./master]
    input = slave
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'plank_right'
    new_block_id = '20'
    new_block_name = 'frictionless_master_subdomain'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    order = ${order}
    block = 'plank block'
    scaling = ${fparse 2.0 / (E_plank + E_block)}
  [../]
  [./disp_y]
    order = ${order}
    block = 'plank block'
    scaling = ${fparse 2.0 / (E_plank + E_block)}
  [../]
  [./temp]
    block = 'plank block'
    initial_condition = 300
    scaling = 1e-5
  [../]
  [./thermal_lm]
    block = 'frictionless_slave_subdomain'
    scaling = 1e-5
  [../]
[]

[Kernels]
  [./diff]
    type = HeatConduction
    variable = temp
    block = 'plank block'
  [../]
  [./dt]
    type = HeatConductionTimeDerivative
    variable = temp
    density_name = 1000
    block = 'plank block'
  [../]
  # [./source]
  #   type = BodyForce
  #   variable = temp
  #   block = 'plank'
  #   value = 1e6
  # [../]
[]

[Modules/TensorMechanics/Master]
  [./action]
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress hydrostatic_stress strain_xx strain_yy strain_zz'
    block = 'plank block'
  [../]
[]

[Constraints]
  [./ced]
    type = GapConductanceConstraint
    variable = thermal_lm
    slave_variable = temp
    k = 10
    use_displaced_mesh = true
    master_boundary = plank_right
    slave_boundary = block_left
    master_subdomain = frictionless_master_subdomain
    slave_subdomain = frictionless_slave_subdomain
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = plank_left
    value = 0.0
  [../]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    boundary = plank_bottom
    value = 0.0
  [../]
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = block_right
    function = '-0.1*sin(4*(t+1.5))+0.02+0.05'
  [../]
  [./right_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = block_right
    function = '-t'
  [../]
  [./right_temp]
    type = DirichletBC
    variable = temp
    boundary = block_right
    value = 300
  [../]
  [./left_temp]
    type = DirichletBC
    variable = temp
    boundary = 'plank_left'
    value = 400
  [../]
[]

[Materials]
  [./plank]
    type = ComputeIsotropicElasticityTensor
    block = 'plank'
    poissons_ratio = 0.3
    youngs_modulus = ${E_plank}
  [../]
  [./block]
    type = ComputeIsotropicElasticityTensor
    block = 'block'
    poissons_ratio = 0.3
    youngs_modulus = ${E_block}
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 'plank block'
  [../]
  [./block_temp]
    type = HeatConductionMaterial
    thermal_conductivity = 1000
    specific_heat = 5
    block = 'block'
  [../]
  [./plank_temp]
    type = HeatConductionMaterial
    thermal_conductivity = 2000
    specific_heat = 10
    block = 'plank'
  [../]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       1e-5          NONZERO               1e-15'
  end_time = 13.5
  dt = 0.1
  dtmin = 0.1
  timestep_tolerance = 1e-6
  line_search = 'contact'
  # automatic_scaling = true
  # compute_scaling_once = false
[]

[Postprocessors]
  [./nl_its]
    type = NumNonlinearIterations
  [../]
  [./total_nl_its]
    type = CumulativeValuePostprocessor
    postprocessor = nl_its
  [../]
  [./l_its]
    type = NumLinearIterations
  [../]
  [./total_l_its]
    type = CumulativeValuePostprocessor
    postprocessor = l_its
  [../]
  [./avg_hydro]
    type = ElementAverageValue
    variable = hydrostatic_stress
    block = 'block'
  [../]
  [./max_hydro]
    type = ElementExtremeValue
    variable = hydrostatic_stress
    block = 'block'
  [../]
  [./min_hydro]
    type = ElementExtremeValue
    variable = hydrostatic_stress
    block = 'block'
    value_type = min
  [../]
  [./avg_vonmises]
    type = ElementAverageValue
    variable = vonmises_stress
    block = 'block'
  [../]
  [./max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    block = 'block'
  [../]
  [./min_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    block = 'block'
    value_type = min
  [../]
  [./block_temp_right_avg]
    type = SideAverageValue
    variable = temp
    boundary = 'block_right'
  [../]
  [./block_temp_left_avg]
    type = SideAverageValue
    variable = temp
    boundary = 'block_left'
  [../]
  [./plank_temp_right_avg]
    type = SideAverageValue
    variable = temp
    boundary = 'plank_right'
  [../]
  [./plank_temp_left_avg]
    type = SideAverageValue
    variable = temp
    boundary = 'plank_left'
  [../]
[]

[Outputs]
  exodus = true
  file_base = ${name}
  # [./comp]
  #   type = CSV
  #   show = 'contact'
  # [../]
  [./out]
    type = CSV
    file_base = '${name}_out'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]
