function replicate_scott,l1b,albedo=albedo,ssa=ssa,sva=sva,sza=sza
  albedo=level_1b_compact(l1b,'albedo')
  ssa=level_1b_compact(l1b,'scattering_angle')
  sva=level_1b_compact(l1b,'view_angle')
  sza=level_1b_compact(l1b,'zenith_angle')
end
