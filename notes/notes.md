
- nakreslit to analyticke riesenie
- a porovnat s numerickym (nejak smart)
- zohnat to analyticke riesenie
- skontrolovat assembly advection termu

TODO
- solutions naming convention
- pridat vypisovacky aj na teplotu
- nakreslit rozne riesenia na porovnanie
- spomenut explicitne v readme ake rovnice to riesi

otazky
- material_data? it?
- preco local_rhs sa pocita pri assemble_system i assemble_rhs? a preco je tak zlozite?

teoria
- explain (flux continuity is enforced weakly by the standard diffusion bilinear form) - heat flux from context
- (steady) advection - diffusion equations - how does the heat transfer work in moving fluid

refactor
- ConjugateHeatTransferSolver is a big class - does it make sense to spolit somehow?
- split main, run to smaller functions (prints, loads)
- no-slip just special case of dirichlet but treated independently and moreover differently, less generally (setup_flow_dofs)

benchmarking
- manufactured solutions
- rozsirit funkcionalitu, aby sa to dalo (source term)

fem-bookmark 2.2.6