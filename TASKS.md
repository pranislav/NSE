1. Keep the current `hp` flow formulation for velocity-pressure and the existing fluid-solid interface no-slip constraints.
2. Add a separate scalar temperature problem on the whole mesh with its own `FE_Q`, `DoFHandler`, constraints, matrix, right-hand side, and solution vectors.
3. Implement phase one as steady diffusion only:
   use a continuous temperature field over both fluid and solid;
   select thermal conductivity by `cell->material_id()`;
   do not add explicit interior interface terms for perfect thermal contact.
4. Apply only external thermal boundary conditions in phase one.
   For the current code, use Dirichlet conditions on boundaries `10` and `20` with different fixed temperatures and leave the other thermal boundaries natural.
5. Solve and output the temperature field together with the flow field on the current mesh.
6. After phase one is stable, add fluid advection to the temperature equation on fluid cells only, using the converged velocity field as a known coefficient.
7. After that, consider refinement strategy for temperature:
   either keep velocity-based refinement for now;
   or combine flow and temperature error indicators.
8. Leave contact resistance and temperature jumps out of scope for now.
   If those are needed later, add explicit interior-face terms on fluid-solid interfaces.
