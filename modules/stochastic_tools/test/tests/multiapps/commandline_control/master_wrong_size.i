[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Problem]
  type = FEProblem
  solve = false
[]

[Executioner]
  type = Steady
[]

[MultiApps]
  [sub]
    type = FullSolveMultiApp
    positions = '0 0 0
                 1 1 1'
    input_files = 'sub.i'
  []
[]

[Distributions]
  [uniform]
    type = Uniform
    lower_bound = 5
    upper_bound = 10
  []
[]

[Samplers]
  [sample]
    type = MonteCarlo
    num_rows = 10
    distributions = 'uniform'
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[Controls]
  [cmdline]
    type = MultiAppCommandLineControl
    multi_app = sub
    sampler = sample
    arguments = 'Mesh/nx'
  []
[]
