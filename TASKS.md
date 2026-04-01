1. Keep the current `hp` flow formulation for velocity-pressure and the existing fluid-solid interface no-slip constraints.
2. Keep a separate scalar temperature problem on the whole mesh with its own `FE_Q`, `DoFHandler`, constraints, matrix, right-hand side, and solution vectors.
3. Preserve the current conjugate-heat-transfer model structure:
   use one continuous temperature field over both fluid and solid;
   select thermal conductivity by `cell->material_id()`;
   do not add explicit interior interface terms for perfect thermal contact.
4. Document the current thermal boundary conditions explicitly.
   The present setup is intentional:
   `T = 0` on outer wall boundaries `10` and `20`;
   `T = 1` on inflow boundary `60`;
   leave the other thermal boundaries natural.
   This differs from the original phase-one note and is meant to make the impact of advection visible.
5. Rework the temperature integration in small verified steps.
   Step 1: replace the temperature linear solver with one appropriate for the non-symmetric advection-diffusion operator.
   Step 2: remove brittle assumptions in fluid/solid classification and velocity lookup between the two DoFHandlers.
   Step 3: add verification output or checks for thermal boundary conditions and expected temperature bounds.
   Step 4: assess whether advection stabilization is needed and add SUPG or another stabilization if required.
   Step 5: revisit refinement and decide whether to keep velocity-only refinement or combine flow and temperature indicators.
6. Leave contact resistance and temperature jumps out of scope for now.
   If those are needed later, add explicit interior-face terms on fluid-solid interfaces.
