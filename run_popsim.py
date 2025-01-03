from popsim import run_population

key = input('Enter key ("visual" for visual simulation, "run" for plotting only) :')

if key == 'visual':

    run_population(key)

else:

    model,ga,gm = run_population(key)