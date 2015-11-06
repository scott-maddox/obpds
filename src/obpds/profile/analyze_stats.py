import pstats
p = pstats.Stats('stats')

# Print all functions, ordered by name
# p.strip_dirs().sort_stats(-1).print_stats()

# print cumulative runtime of the 10 longest running lines
# p.sort_stats('cumulative').print_stats(20)
 
# print non-cumulative runtime of the 10 longest running lines
p.sort_stats('time').print_stats()