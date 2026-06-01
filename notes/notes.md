- skontrolovat assembly advection termu

TODO
- pridat vypisovacky aj na teplotu
- spomenut explicitne v readme ake rovnice to riesi
- ako spustat generate do README
- compute_errors - zmenit specifikaciu tlakovej komponenty z dim na pressure(enum)
- enum for block indices?
- case files refinement cycles je ich o jeden viac
- implement no-penetration bc option pre vodnu hladinu

otazky
- material_data? it?
- preco local_rhs sa pocita pri assemble_system i assemble_rhs? a preco je tak zlozite?

teoria
- explain (flux continuity is enforced weakly by the standard diffusion bilinear form) - heat flux from context
- (steady) advection - diffusion equations - how does the heat transfer work in moving fluid

refactor
- ConjugateHeatTransferSolver is a big class - does it make sense to spolit somehow?
- split main, run to smaller functions (prints, loads)

benchmarking
- manufactured solutions
