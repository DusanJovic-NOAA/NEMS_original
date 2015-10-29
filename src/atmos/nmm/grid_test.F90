program grid_test
  use module_REDUCTION
  type(grid) :: grid

  grid%ims=1   grid%ids=1
  grid%ime=30  grid%ide=30
  grid%jms=1   grid%jds=1
  grid%jme=40  grid%jde=40
  grid%kms=1   grid%kds=1
  grid%kme=50  grid%kde=50


end program grid_test
