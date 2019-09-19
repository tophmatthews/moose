# Break Mesh By Block

!syntax description /MeshModifiers/BreakMeshByBlock

This class implement a MeshModifiers to split a monolithic mesh by blocks similar to what is proposed by VP Nguyen [!cite](Nguyen2014).

To split the mesh, nodes shared by multiple blocks are duplicated N-1 times (where N is the number of blocks sharing a particular node). Each duplicated nodes is assigned to one block and all the element sharing that node are updated. A new sideset identifying the new interface is added and it is always linked to elements belonging to blocks with the lower id.



As an option, the interface can be split into $Q$ different sidesets. $Q$ is the number of adjacent block pairs. This is achieved by setting  `split_interface=true`. This is useful when modeling interfaces with different parameters.

!syntax inputs /MeshModifiers/BreakMeshByBlock

!syntax parameters /MeshModifiers/BreakMeshByBlock

!bibtex bibliography
