# This test checks the thermal expansion calculated via a mean thermal expansion coefficient.
# The coefficient is selected so as to result in a 1e-4 strain in the x-axis, and to cross over
# from positive to negative strain.

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[AuxVariables]
  [./temp]
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    eigenstrain_names = eigenstrain
    generate_output = 'strain_xx strain_yy strain_zz'
    use_automatic_differentiation = true
  [../]
[]

[BCs]
  [./left]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  [./bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  [./back]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
[]

[AuxKernels]
  [./temp]
    type = FunctionAux
    variable = temp
    function = '1 + t'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ADComputeLinearElasticStress
  [../]
  [./thermal_expansion_strain]
    type = ADComputeMeanThermalExpansionPolynomialEigenstrain
    thermal_expansion_coefficients = '${fparse 1e-4 / (2 - 1.5 - 1e-4 * (1.5 - 1.2))}'
    thermal_expansion_function_reference_temperature = 1.2
    stress_free_temperature = 1.5
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]
[]

[Postprocessors]
  [./disp_x_max]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  [../]
  [./temp_avg]
    type = ElementAverageValue
    variable = temp
  [../]
[]

[Executioner]
  type = Transient

  end_time = 1.0
  dt = 0.1
[]

[Outputs]
  csv = true
[]
