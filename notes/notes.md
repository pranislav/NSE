TODO
- compute_errors - zmenit specifikaciu tlakovej komponenty z dim na pressure(enum) -- cant make enum if dim is not known / mms built for 2D - enum implemented
- enum for block indices?
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
- update teploty inside newton iteration - velmi nepekne
