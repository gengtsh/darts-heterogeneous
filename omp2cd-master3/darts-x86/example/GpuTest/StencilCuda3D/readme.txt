2,2i and 2i.p: cpu without blocking, cutting based on slice
StencilCuda3D_7Points.2:  cuda: row/cols/slice uint64_t
StencilCuda3D_7Points.2i: row/cols/slice int
StencilCuda3D_7Points.2i.p: collect data for for profiler
StencilCuda3D_7Points.2i.s: based on 2i, add IsStatic function

StencilCuda3D_7Points.3i: cpu with blocking, cutting based on any cube; use cuda 3d copy
StencilCuda3D_7Points.4i: based on 3i, use cuda 3d copy
StencilCuda3D_7Points.5i: based on 4i, use cuda 3d copy, add IsStatic function
 