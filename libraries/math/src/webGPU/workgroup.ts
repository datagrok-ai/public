/** Threads per dimension of the 2-D workgroups of this library's compute shaders: 8 x 8 = 64 invocations, a whole
 * number of subgroups on every GPU (subgroups are 8, 16, 32 or 64 wide). A workgroup that is not (10 x 10 = 100, the size
 * until 2026-10-14) ends in a partly filled subgroup, and Intel's Windows driver did not run that subgroup's invocations
 * (96-99 of every workgroup) once the same pipeline was created again in a page: MMP analysis then found no pairs and
 * dispatched forever from its second run on, and any kernel of this library could lose those threads' work. */
export const WORKGROUP_DIM = 8;
