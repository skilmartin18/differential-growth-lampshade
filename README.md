
multiprocessed:


2026-02-06 15:05:26,609 - INFO - total runtime 15.303521637993981
2026-02-06 15:05:26,609 - INFO - avg time per iteration 0.05101179613998587


2026-02-06 15:12:31,502 - INFO - total runtime 48.80206467500102
2026-02-06 15:12:31,502 - INFO - avg time per iteration 0.0976041617699957

1000 its
2026-02-06 15:35:33,052 - INFO - total runtime 358.5061117239529
2026-02-06 15:35:33,052 - INFO - avg time per iteration 0.3585061117239529

Non-multiprocessed:
300 its
2026-02-06 15:06:58,065 - INFO - total runtime 19.842740174004575
2026-02-06 15:06:58,066 - INFO - avg time per iteration 0.06614252107334323

500 its
2026-02-06 15:10:05,040 - INFO - total runtime 106.61456236900267
2026-02-06 15:10:05,040 - INFO - avg time per iteration 0.21322915168800682





- make asymmetry work better DONE
- make random perturbations not break the sim DONE
- implement cell death DONE works well
- implement random cell life-times
- implement cell enhanced growing random/preset
- implement growing regions -> position in space based (depending on region of space the cell inhabits)
- implement growing regions -> pick a group of consecutive cells and increase their lifespan, or friction, or growth multiplier? preset/random
- implement sun (single point of attractive force)
- implement moon (rotating point of attractive force)


- simulation -> make save files for runs, periodically save logs and image of cells 